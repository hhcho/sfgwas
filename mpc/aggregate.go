package mpc

import (
	// "bufio"

	"github.com/hhcho/sfgwas/crypto"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"

	"github.com/ldsec/lattigo/v2/dckks"
)

func (netObj *Network) AggregateDecryptSharesMat(poly [][]dckks.PCKSShare, outLevel int) (out [][]dckks.PCKSShare) {
	out = make([][]dckks.PCKSShare, len(poly))
	pid := netObj.pid
	for i := range out {
		out[i] = make([]dckks.PCKSShare, len(poly[i]))
	}

	if pid > 0 {
		if pid == netObj.hubPid {
			for r := range poly {
				for c := range poly[r] {
					for i := range poly[r][c] {
						out[r][c][i] = netObj.dckksContext.RingQ.NewPolyLvl(outLevel)
					}
				}
			}
			for p := 1; p < netObj.NumParties; p++ {
				for r := range poly {
					for c := range poly[r] {
						for i := range poly[r][c] {
							var newPoly *ring.Poly
							if p != pid {
								newPoly = netObj.ReceivePoly(p)
							} else {
								newPoly = poly[r][c][i]
							}
							level := len(newPoly.Coeffs) - 1
							netObj.dckksContext.RingQ.AddLvl(level, newPoly, out[r][c][i], out[r][c][i])
						}
					}
				}
			}
			for p := 1; p < netObj.NumParties; p++ {
				if p != pid {
					for r := range poly {
						for c := range poly[r] {
							for i := range poly[r][c] {
								netObj.SendPoly(out[r][c][i], p)
							}
						}
					}
				}
			}
		} else {
			for r := range poly {
				for c := range poly[r] {
					for i := range poly[r][c] {
						netObj.SendPoly(poly[r][c][i], netObj.hubPid)
					}
				}
			}
			for r := range poly {
				for c := range poly[r] {
					for i := range poly[r][c] {
						out[r][c][i] = netObj.ReceivePoly(netObj.hubPid)
					}
				}
			}
		}
	}
	return
}

// level := uint64(len(share1[0].Coeffs)) - 1
// pcks.dckksContext.RingQ.AddLvl(level, share1[0], share2[0], shareOut[0])
// pcks.dckksContext.RingQ.AddLvl(level, share1[1], share2[1], shareOut[1])
func (netObj *Network) AggregateDecryptShares(poly *dckks.PCKSShare, outLevel int) (out *dckks.PCKSShare) {
	out = new(dckks.PCKSShare)

	pid := netObj.pid
	if pid > 0 {
		if pid == netObj.hubPid {
			for i := range poly {
				out[i] = netObj.dckksContext.RingQ.NewPolyLvl(outLevel)
			}
			for p := 1; p < netObj.NumParties; p++ {
				for i := range poly {
					var newPoly *ring.Poly
					if p != pid {
						newPoly = netObj.ReceivePoly(p)
					} else {
						newPoly = poly[i]
					}
					level := len(newPoly.Coeffs) - 1
					netObj.dckksContext.RingQ.AddLvl(level, newPoly, out[i], out[i])
				}
			}
			for p := 1; p < netObj.NumParties; p++ {
				if p != pid {
					for i := range poly {
						netObj.SendPoly(out[i], p)
					}
				}
			}
		} else {
			for i := range poly {
				netObj.SendPoly(poly[i], netObj.hubPid)
			}
			for i := range poly {
				out[i] = netObj.ReceivePoly(netObj.hubPid)
			}
		}
	}
	return
}

func (netObj *Network) AggregatePubKeyShares(poly *drlwe.CKGShare) (out *drlwe.CKGShare) {
	out = new(drlwe.CKGShare)

	pid := netObj.pid
	if pid > 0 {
		if pid == netObj.hubPid {
			out.Poly = netObj.dckksContext.RingQP.NewPoly()
			for p := 1; p < netObj.NumParties; p++ {
				var newPoly *ring.Poly
				if p != pid {
					newPoly = netObj.ReceivePoly(p)
				} else {
					newPoly = poly.Poly
				}
				netObj.dckksContext.RingQP.Add(newPoly, out.Poly, out.Poly)
			}
			for p := 1; p < netObj.NumParties; p++ {
				if p != pid {
					netObj.SendPoly(out.Poly, p)
				}
			}
		} else {
			netObj.SendPoly(poly.Poly, netObj.hubPid)
			out.Poly = netObj.ReceivePoly(netObj.hubPid)
		}
	}
	return
}

func (netObj *Network) AggregateRelinKeyShare(share *drlwe.RKGShare, secondSlot bool) (shareOut *drlwe.RKGShare) {
	dckksContext := netObj.dckksContext
	pid := netObj.pid
	contextQP := dckksContext.RingQP

	shareOut = new(drlwe.RKGShare)
	shareOut.SetValue(make([][2]*ring.Poly, dckksContext.Beta))

	if pid > 0 {
		if pid == netObj.hubPid {
			// Initialize
			for i := 0; i < dckksContext.Beta; i++ {
				shareOut.Value()[i][0] = contextQP.NewPoly()
				if secondSlot {
					shareOut.Value()[i][1] = contextQP.NewPoly()
				}
			}

			// Aggregate
			for p := 1; p < netObj.NumParties; p++ {
				for i := 0; i < dckksContext.Beta; i++ {
					var other0, other1 *ring.Poly
					if p != pid {
						other0 = netObj.ReceivePoly(p)
						if secondSlot {
							other1 = netObj.ReceivePoly(p)
						}
					} else {
						other0 = share.Value()[i][0]
						if secondSlot {
							other1 = share.Value()[i][1]
						}
					}
					contextQP.Add(other0, shareOut.Value()[i][0], shareOut.Value()[i][0])
					if secondSlot {
						contextQP.Add(other1, shareOut.Value()[i][1], shareOut.Value()[i][1])
					}
				}
			}

			// Broadcast
			for p := 1; p < netObj.NumParties; p++ {
				if p != pid {
					for i := 0; i < dckksContext.Beta; i++ {
						netObj.SendPoly(shareOut.Value()[i][0], p)
						if secondSlot {
							netObj.SendPoly(shareOut.Value()[i][1], p)
						}
					}
				}
			}
		} else {
			// Send share
			for i := 0; i < dckksContext.Beta; i++ {
				netObj.SendPoly(share.Value()[i][0], netObj.hubPid)
				if secondSlot {
					netObj.SendPoly(share.Value()[i][1], netObj.hubPid)
				}
			}

			// Receive result
			for i := 0; i < dckksContext.Beta; i++ {
				shareOut.Value()[i][0] = netObj.ReceivePoly(netObj.hubPid)
				if secondSlot {
					shareOut.Value()[i][1] = netObj.ReceivePoly(netObj.hubPid)
				}
			}
		}
	}

	return
}

func (netObj *Network) AggregateRotKeyShare(share *drlwe.RTGShare) (shareOut *drlwe.RTGShare) {
	inShare := new(drlwe.RKGShare)
	inShare.SetValue(make([][2]*ring.Poly, len(share.Value)))
	for i := range inShare.Value() {
		inShare.Value()[i][0] = share.Value[i]
	}

	out := netObj.AggregateRelinKeyShare(inShare, false).Value()

	shareOut = new(drlwe.RTGShare)
	shareOut.Value = make([]*ring.Poly, len(out))
	for i := range out {
		shareOut.Value[i] = out[i][0]
	}

	return
}

func (netObj *Network) AggregateRefreshShareMat(share [][]ring.Poly, outLevel int) (shareOut [][]ring.Poly) {
	contextQ := netObj.dckksContext.RingQ
	shareOut = make([][]ring.Poly, len(share))
	pid := netObj.pid
	if pid > 0 {
		if pid == netObj.hubPid {
			// Initialize
			for i := range shareOut {
				shareOut[i] = make([]ring.Poly, len(share[i]))
				for j := range shareOut[i] {
					shareOut[i][j] = *contextQ.NewPolyLvl(outLevel)
				}

			}

			// Aggregate
			for p := 1; p < netObj.NumParties; p++ {
				var other [][]ring.Poly
				//send all at once
				if p != pid {
					other = netObj.ReceivePolyMat(p)
				} else {
					other = share
				}

				for i := range share {
					for j := range share[i] {
						contextQ.AddLvl(len(other[i][j].Coeffs)-1, &other[i][j], &shareOut[i][j], &shareOut[i][j])
					}

				}
			}

			// Broadcast
			for p := 1; p < netObj.NumParties; p++ {
				if p != pid {
					netObj.SendPolyMat(shareOut, p)
				}
			}
		} else {
			// Send share
			netObj.SendPolyMat(share, netObj.hubPid)
			shareOut = netObj.ReceivePolyMat(netObj.hubPid)

		}
	}
	return

}

func (netObj *Network) AggregateRefreshShareVec(share []*ring.Poly, outLevel int) (shareOut []*ring.Poly) {
	contextQ := netObj.dckksContext.RingQ
	shareOut = make([]*ring.Poly, len(share))
	pid := netObj.pid
	if pid > 0 {
		if pid == netObj.hubPid {
			// Initialize
			for i := range shareOut {
				shareOut[i] = contextQ.NewPolyLvl(outLevel)
			}

			// Aggregate
			for p := 1; p < netObj.NumParties; p++ {
				var other *ring.Poly
				for i := range share {
					if p != pid {
						other = netObj.ReceivePoly(p)
					} else {
						other = share[i]
					}
					contextQ.AddLvl(len(other.Coeffs)-1, other, shareOut[i], shareOut[i])
				}
			}

			// Broadcast
			for p := 1; p < netObj.NumParties; p++ {
				if p != pid {
					for i := range share {
						netObj.SendPoly(shareOut[i], p)
					}
				}
			}
		} else {
			// Send share
			for i := range share {
				netObj.SendPoly(share[i], netObj.hubPid)
			}

			// Receive result
			for i := range share {
				shareOut[i] = netObj.ReceivePoly(netObj.hubPid)
			}
		}
	}
	return
}

func (netObj *Network) AggregateRefreshShare(share *ring.Poly, outLevel int) (shareOut *ring.Poly) {
	contextQ := netObj.dckksContext.RingQ
	pid := netObj.pid
	if pid > 0 {
		if pid == netObj.hubPid {
			// Initialize
			shareOut = contextQ.NewPolyLvl(outLevel)

			// Aggregate
			for p := 1; p < netObj.NumParties; p++ {
				var other *ring.Poly
				if p != pid {
					other = netObj.ReceivePoly(p)
				} else {
					other = share
				}
				contextQ.AddLvl(len(other.Coeffs)-1, other, shareOut, shareOut)
			}

			// Broadcast
			for p := 1; p < netObj.NumParties; p++ {
				if p != pid {
					netObj.SendPoly(shareOut, p)
				}
			}
		} else {
			// Send share
			netObj.SendPoly(share, netObj.hubPid)

			// Receive result
			shareOut = netObj.ReceivePoly(netObj.hubPid)
		}
	}
	return
}

func (netObj *Network) AggregateCText(cryptoParams *crypto.CryptoParams, val *ckks.Ciphertext) (out *ckks.Ciphertext) {
	pid := netObj.GetPid()
	if pid == 0 {
		return nil
	}

	if pid == netObj.GetHubPid() {
		//receive and add
		out = val.CopyNew().Ciphertext()
		cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			for p := 1; p < netObj.GetNParty(); p++ {
				if p != pid {
					other := netObj.ReceiveCiphertext(cryptoParams, p)
					eval.Add(other, out, out)
				}
			}
			return nil
		})

		for p := 1; p < netObj.GetNParty(); p++ {
			if p != pid {
				netObj.SendCiphertext(out, p)
			}
		}
	} else {
		netObj.SendCiphertext(val, netObj.GetHubPid())
		out = netObj.ReceiveCiphertext(cryptoParams, netObj.GetHubPid())
	}

	return
}

func (netObj *Network) AggregateIntVec(vec []uint64) (out []uint64) {
	pid := netObj.GetPid()
	if pid == 0 {
		return nil
	}

	if pid == netObj.GetHubPid() {
		//receive and add
		out = make([]uint64, len(vec))
		copy(out, vec)

		for p := 1; p < netObj.GetNParty(); p++ {
			if p != pid {
				other := netObj.ReceiveIntVector(len(vec), p)
				for i := range other {
					out[i] += other[i]
				}
			}
		}

		for p := 1; p < netObj.GetNParty(); p++ {
			if p != pid {
				netObj.SendIntVector(out, p)
			}
		}
	} else {
		netObj.SendIntVector(vec, netObj.GetHubPid())
		out = netObj.ReceiveIntVector(len(vec), netObj.GetHubPid())
	}

	return
}

func (netObj *Network) AggregateSharesCT(cryptoParams *crypto.CryptoParams, ct *ckks.Ciphertext) *ckks.Ciphertext {
	agg := ct.CopyNew().Ciphertext()

	//get shares from everyone except for pid = 0
	pid := netObj.GetPid()
	for i := 1; i < netObj.GetNParty(); i++ {
		if pid == i {
			continue
		} else if i > pid {
			//receive first and then send
			agg = crypto.Add(cryptoParams, agg, netObj.ReceiveCiphertext(cryptoParams, i))
			netObj.SendCiphertext(ct, i)
		} else {
			//send and then receive
			netObj.SendCiphertext(ct, i)
			agg = crypto.Add(cryptoParams, agg, netObj.ReceiveCiphertext(cryptoParams, i))
		}

	}
	return agg

}

func (netObj *Network) AggregateCVec(cryptoParams *crypto.CryptoParams, vec crypto.CipherVector) crypto.CipherVector {
	return netObj.AggregateCMat(cryptoParams, crypto.CipherMatrix{vec})[0]
}

func (netObj *Network) AggregateCMat(cryptoParams *crypto.CryptoParams, mat crypto.CipherMatrix) (out crypto.CipherMatrix) {
	pid := netObj.GetPid()
	if pid == 0 {
		return nil
	}

	if pid == netObj.GetHubPid() {
		//receive and add
		out = crypto.CopyEncryptedMatrix(mat)
		cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			for p := 1; p < netObj.GetNParty(); p++ {
				if p != pid {
					other := netObj.ReceiveCipherMatrix(cryptoParams, len(out), len(out[0]), p)
					for i := range other {
						for j := range other[i] {
							eval.Add(other[i][j], out[i][j], out[i][j])
						}
					}
				}
			}
			return nil
		})

		for p := 1; p < netObj.GetNParty(); p++ {
			if p != pid {
				netObj.SendCipherMatrix(out, p)
			}
		}
	} else {
		netObj.SendCipherMatrix(mat, netObj.GetHubPid())
		out = netObj.ReceiveCipherMatrix(cryptoParams, len(mat), len(mat[0]), netObj.GetHubPid())
	}

	return
}
