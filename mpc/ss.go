package mpc

import (
	"math/big"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/ldsec/lattigo/v2/utils"

	// "go.dedis.ch/onet/v3/log"

	"github.com/hhcho/sfgwas-private/crypto"

	"github.com/ldsec/lattigo/v2/dckks"
	"github.com/ldsec/lattigo/v2/ring"

	//"github.com/hhcho/mpc-core"
	"github.com/ldsec/lattigo/v2/ckks"
	//"math/bits"
)

func (mpcObj *MPC) SSMultMat(a, b mpc_core.RMat) mpc_core.RMat {
	ar, am := mpcObj.BeaverPartitionMat(a)
	br, bm := mpcObj.BeaverPartitionMat(b)
	ab := mpcObj.BeaverMultMat(ar, am, br, bm)
	return mpcObj.BeaverReconstructMat(ab)
}

func (mpcObj *MPC) SSMultElemVecScalar(a mpc_core.RVec, b mpc_core.RElem) mpc_core.RVec {
	ar, am := mpcObj.BeaverPartitionVec(a)
	br, bm := mpcObj.BeaverPartition(b)
	x := mpc_core.InitRVec(mpcObj.rtype.Zero(), len(a))
	for i := range x {
		x[i] = mpcObj.BeaverMult(ar[i], am[i], br, bm)
	}
	return mpcObj.BeaverReconstructVec(x)
}

func (mpcObj *MPC) SSSquareElemVec(a mpc_core.RVec) mpc_core.RVec {
	ar, am := mpcObj.BeaverPartitionVec(a)
	x := mpcObj.BeaverMultElemVec(ar, am, ar, am)
	return mpcObj.BeaverReconstructVec(x)
}

func (mpcObj *MPC) SSMultElemVec(a, b mpc_core.RVec) mpc_core.RVec {
	ar, am := mpcObj.BeaverPartitionVec(a)
	br, bm := mpcObj.BeaverPartitionVec(b)
	x := mpcObj.BeaverMultElemVec(ar, am, br, bm)
	return mpcObj.BeaverReconstructVec(x)
}

func (mpcObj *MPC) SSMultElemMat(a, b mpc_core.RMat) mpc_core.RMat {
	ar, am := mpcObj.BeaverPartitionMat(a)
	br, bm := mpcObj.BeaverPartitionMat(b)
	x := mpcObj.BeaverMultElemMat(ar, am, br, bm)
	return mpcObj.BeaverReconstructMat(x)
}

// Row-major
func (mpcObj *MPC) SSToCMat(cryptoParams *crypto.CryptoParams, rm mpc_core.RMat) (cm crypto.CipherMatrix) {
	if mpcObj.GetPid() == 0 {
		cm = make(crypto.CipherMatrix, 1)
		cm[0] = make(crypto.CipherVector, 1)
		return
	}

	rtype := rm.Type().Zero()

	if rtype.TypeID() != mpc_core.LElem256UniqueID && rtype.TypeID() != mpc_core.LElem128UniqueID {
		panic("SSToCMat only supported for LElem128 or LElem256")
	}

	slots := cryptoParams.GetSlots()
	numCtxRow := len(rm)
	nElemCol := len(rm[0])
	numCtxCol := 1 + ((nElemCol - 1) / slots)

	// Pad RVec
	//rvNew := mpc_core.InitRVec(rtype.Zero(), cryptoParams.GetSlots()*numCtxCol)
	//for i := range rvNew {
	//	rvNew[i] = rv[i % len(rv)]
	//}
	//rv = rvNew

	bound := rtype.Modulus()
	bound.Quo(bound, big.NewInt(4*int64(mpcObj.GetNParty()-1)))
	boundHalf := new(big.Int).Rsh(bound, 1)

	mask := make(mpc_core.RMat, len(rm))
	boundElem := rtype.FromBigInt(bound)
	for i := range rm {
		mask[i] = make(mpc_core.RVec, len(rm[0]))
		for j := range rm[0] {
			// log.LLvl1(time.Now().Format(time.RFC3339), "ss to cvec bound: ", bound)
			tmp := ring.RandInt(bound)
			mask[i][j] = rtype.FromBigInt(tmp)
			if tmp.Cmp(boundHalf) >= 0 {
				mask[i][j] = mask[i][j].Sub(boundElem)
			}
		}
	}

	rmMask := rm.Copy()
	rmMask.Sub(mask)
	rmMask = mpcObj.RevealSymMat(rmMask)

	var share mpc_core.RMat
	if mpcObj.GetPid() == mpcObj.GetHubPid() { // share = (x - r) + r_i
		share = rmMask
		share.Add(mask)
	} else { // share = r_i
		share = mask
	}

	pm := make(crypto.PlainMatrix, numCtxRow)
	for i := range pm {
		pm[i] = make(crypto.PlainVector, numCtxCol)
		cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
			start := 0
			end := slots
			for j := 0; j < numCtxCol; j++ {
				if end > nElemCol {
					end = nElemCol
				}

				pm[i][j] = encoder.EncodeRVecNew(share[i][start:end], uint64(end-start), mpcObj.GetFracBits())

				start += slots
				end += slots
			}
			return nil
		})
	}

	cm = crypto.EncryptPlaintextMatrix(cryptoParams, pm)

	return mpcObj.Network.AggregateCMat(cryptoParams, cm)
}

func (mpcObj *MPC) SSToCVec(cryptoParams *crypto.CryptoParams, rv mpc_core.RVec) (cv crypto.CipherVector) {
	return mpcObj.SSToCMat(cryptoParams, mpc_core.RMat{rv})[0]
}
func (mpcObj *MPC) SStoCiphertext(cryptoParams *crypto.CryptoParams, rv mpc_core.RVec) *ckks.Ciphertext {
	return mpcObj.SSToCVec(cryptoParams, rv)[0]
}

func (mpcObj *MPC) CMatToSS(cryptoParams *crypto.CryptoParams, rtype mpc_core.RElem, cm crypto.CipherMatrix, sourcePid, numCtxRow, numCtxCol, nElemRow int) (rm mpc_core.RMat) {
	slots := cryptoParams.GetSlots()
	fracBits := mpcObj.GetFracBits()

	rm = mpc_core.InitRMat(rtype.Zero(), numCtxRow, nElemRow)
	if mpcObj.GetPid() == 0 {
		return
	}

	if sourcePid > 0 {
		cm = mpcObj.Network.BroadcastCMat(cryptoParams, cm, sourcePid, numCtxRow, numCtxCol)
	}
	cm, levelStart := crypto.FlattenLevels(cryptoParams, cm)
	ctMask := crypto.CopyEncryptedMatrix(cm)

	paramN := cryptoParams.Params.N()

	dckksContext := dckks.NewContext(cryptoParams.Params)
	shareDecrypt := make([][]*ring.Poly, numCtxRow)
	for i := range shareDecrypt {
		shareDecrypt[i] = make([]*ring.Poly, numCtxCol)
		for j := range shareDecrypt[i] {
			shareDecrypt[i][j] = dckksContext.RingQ.NewPolyLvl(levelStart)
		}
	}
	maskBigint := make([][][]*big.Int, numCtxRow)
	for i := range maskBigint {
		maskBigint[i] = make([][]*big.Int, numCtxCol)
		for j := range maskBigint[i] {
			maskBigint[i][j] = make([]*big.Int, paramN)
		}
	}

	context := dckksContext.RingQ

	prng, _ := utils.NewPRNG()
	sampler := ring.NewGaussianSampler(prng)

	bound := ring.NewUint(context.Modulus[0])
	for i := 1; i < levelStart+1; i++ {
		bound.Mul(bound, ring.NewUint(context.Modulus[i]))
	}
	bound.Quo(bound, ring.NewUint(2*uint64(mpcObj.GetNParty()-1)))

	//fmt.Println("CMatToSS: Bound bit length ", bound.BitLen())

	// Check if there is enough space in ct for masks
	//if bound.Cmp(##) < 0 {
	//	panic(fmt.Sprintf("Attempted SS conversion on a ciphertext without enough levels -> %d", levelStart))
	//}

	boundHalf := new(big.Int).Rsh(bound, 1)

	var sign int
	for k := range cm {
		for i := range cm[k] {
			for j := range maskBigint[k][i] {
				// TODO: check relation between coeff size and decoded output size
				m := ring.RandInt(bound)
				sign = m.Cmp(boundHalf)
				if sign == 1 || sign == 0 {
					m.Sub(m, bound)
				}
				maskBigint[k][i][j] = m
			}
		}
	}

	for k := range cm {
		for i := range cm[k] {
			// h0 = mask (at level min)
			context.SetCoefficientsBigintLvl(levelStart, maskBigint[k][i], shareDecrypt[k][i])
			context.NTTLvl(levelStart, shareDecrypt[k][i], shareDecrypt[k][i])

			ctMask[k][i].SetValue([]*ring.Poly{shareDecrypt[k][i].CopyNew(), context.NewPoly()})

			// h0 = sk*c1 + mask
			context.MulCoeffsMontgomeryAndAddLvl(levelStart, cryptoParams.Sk.Value, cm[k][i].Value()[1], shareDecrypt[k][i])

			// h0 = sk*c1 + mask + e0
			tmp := sampler.ReadNew(dckksContext.RingQ, 3.19, 19)
			dckksContext.RingQ.NTT(tmp, tmp)

			context.AddLvl(levelStart, shareDecrypt[k][i], tmp, shareDecrypt[k][i])
		}
	}

	// TODO: communicate in one batch
	agg := make([][]*ring.Poly, len(shareDecrypt))
	for i := range shareDecrypt {
		agg[i] = mpcObj.Network.AggregateRefreshShareVec(shareDecrypt[i], levelStart)
	}

	pt := make(crypto.PlainMatrix, numCtxRow)
	ptMask := make(crypto.PlainMatrix, numCtxRow)
	for i := range cm {
		pt[i] = make(crypto.PlainVector, numCtxCol)
		ptMask[i] = make(crypto.PlainVector, numCtxCol)

		for j := range cm[i] {
			ctOut := cm[i][j].CopyNew().Ciphertext()
			context.AddLvl(levelStart, ctOut.Value()[0], agg[i][j], ctOut.Value()[0])

			pt[i][j] = ctOut.Plaintext()
			ptMask[i][j] = ctMask[i][j].Plaintext()
		}
	}

	rm = mpc_core.InitRMat(rtype.Zero(), numCtxRow, nElemRow)
	cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
		for i := range pt {
			for j := range pt[i] {
				var rvOut mpc_core.RVec
				if mpcObj.GetPid() == mpcObj.GetHubPid() {
					rvOut = encoder.DecodeRVec(rtype, pt[i][j], uint64(slots), fracBits)
				} else {
					rvOut = mpc_core.InitRVec(rtype.Zero(), slots)
				}
				rvMask := encoder.DecodeRVec(rtype, ptMask[i][j], uint64(slots), fracBits)
				rvOut.Sub(rvMask)

				start := j * slots
				end := start + slots
				if end > nElemRow {
					end = nElemRow
				}

				for k := 0; k < end-start; k++ {
					rm[i][start+k] = rvOut[k]
				}
			}
		}
		return nil
	})
	return
}

func (mpcObj *MPC) CVecToSS(cryptoParams *crypto.CryptoParams, rtype mpc_core.RElem, cv crypto.CipherVector, sourcePid, numCtx, nElem int) (rm mpc_core.RVec) {
	return mpcObj.CMatToSS(cryptoParams, rtype, crypto.CipherMatrix{cv}, sourcePid, 1, numCtx, nElem)[0]
}

func (mpcObj *MPC) CiphertextToSS(cryptoParams *crypto.CryptoParams, rtype mpc_core.RElem, ct *ckks.Ciphertext, sourcePid, N int) (rv mpc_core.RVec) {
	return mpcObj.CVecToSS(cryptoParams, rtype, crypto.CipherVector{ct}, sourcePid, 1, N)
}
