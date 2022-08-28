package gwas

import (
	"fmt"
	"math"
	"time"

	"go.dedis.ch/onet/v3/log"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas-private/crypto"

	"github.com/ldsec/lattigo/v2/ckks"
)

type PCA struct {
	general *ProtocolInfo

	// Data files for geno input
	genoRed  *GenoFileStream
	genoRedT *GenoFileStream

	// Cache files for geno input (to be created)
	genoRedCache  string
	genoRedTCache string

	numInds []int
	numSnps int

	Q crypto.CipherMatrix
}

func (g *ProtocolInfo) InitPCA(genoRed, genoRedT *GenoFileStream) *PCA {

	return &PCA{
		general:  g,
		genoRed:  genoRed,
		genoRedT: genoRedT,

		genoRedCache:  g.CachePath("geno_pca_diag"),
		genoRedTCache: g.CachePath("geno_pca_transpose_diag"),

		numInds: g.gwasParams.FiltNumInds(),
		numSnps: g.gwasParams.numSnpsPCA,
	}

}

func (pca *PCA) DistributedPCA() crypto.CipherMatrix {
	debug := pca.general.config.Debug
	restartIter := pca.general.config.PCARestartIter
	skipPowerIter := pca.general.config.SkipPowerIter
	binaryVersion := pca.general.config.MpcBooleanShares

	gwasParams := pca.general.gwasParams

	X := pca.genoRed   // individual by SNP
	XT := pca.genoRedT // SNP by individual

	Xcache := pca.genoRedCache   // individual by SNP
	XTcache := pca.genoRedTCache // SNP by individual

	nRowsAll := pca.numInds
	cryptoParams := pca.general.cps

	mpcPar := pca.general.mpcObj
	mpcObj := pca.general.mpcObj[0]
	pid := mpcObj.GetPid()

	nsnp, nind := pca.numSnps, 0
	if pid > 0 {
		nind = int(X.NumRows())
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Distributed PCA called: nsnp", nsnp, "nind", nind)

	rtype := mpcObj.GetRType().Zero()
	fracBits := mpcObj.GetFracBits()
	dataBits := mpcObj.GetDataBits()
	slots := cryptoParams.GetSlots()
	totInd := 0
	for i := range nRowsAll {
		totInd += nRowsAll[i]
	}

	npc := gwasParams.NumPC()
	kp := gwasParams.NumPC() + pca.general.config.NumOversample
	nPowerIter := pca.general.config.NumPowerIters

	numSnpSqrtInv := 1.0 / math.Sqrt(float64(nsnp))
	numTotIndSqrtInv := 1.0 / math.Sqrt(float64(totInd))

	// Mean, stdev calculation
	xsum := make([]uint64, nsnp)
	x2sum := make([]uint64, nsnp)
	bucketCount := make([]uint64, kp)
	posCount := make([]uint64, kp)

	// Sketching (nsnp x nind --> nsnp x kp)
	localSketch := make([][]float64, kp)
	for i := range localSketch {
		localSketch[i] = make([]float64, nsnp)
	}

	Zmat := mpc_core.InitRMat(rtype.Zero(), kp, kp)
	var Q, Qloc crypto.CipherMatrix
	var XMean, XVar, XStdInv crypto.CipherVector

	// Preprocess X
	if pid > 0 {
		log.LLvl1(time.Now().Format(time.RFC3339), "Preprocessing X")
		MatMult4StreamPreprocess(cryptoParams, X, 5, Xcache)
		MatMult4StreamPreprocess(cryptoParams, XT, 5, XTcache)
	}

	sx := mpc_core.InitRVec(rtype.Zero(), nsnp)
	sx2 := mpc_core.InitRVec(rtype.Zero(), nsnp)

	log.LLvl1(time.Now().Format(time.RFC3339), "Before sketch")
	mpcObj.AssertSync()

	if pid > 0 {
		log.LLvl1(time.Now().Format(time.RFC3339), "Sketching")

		randIndex := make([]int, nind)
		sgn := make([]float64, nind)

		for i := range randIndex {

			randIndex[i] = mpcObj.Network.Rand.CurPRG().Intn(kp)
			sgn[i] = float64(mpcObj.Network.Rand.CurPRG().Intn(2)*2 - 1)

			bucketCount[randIndex[i]]++
			if sgn[i] > 0 {
				posCount[randIndex[i]]++
			}
		}

		if debug {
			sgnInt := make([]int, nind)
			for i := range sgnInt {
				if sgn[i] > 0 {
					sgnInt[i] = 1
				} else {
					sgnInt[i] = -1
				}
			}
			SaveIntVectorToFile(pca.general.CachePath("SketchSign.txt"), sgnInt)
			SaveIntVectorToFile(pca.general.CachePath("SketchBucketId.txt"), randIndex)
		}

		X.Reset()
		for i := 0; i < nind; i++ {
			row := X.NextRow()

			for j := range row {
				localSketch[randIndex[i]][j] += sgn[i] * float64(row[j])

				xsum[j] += uint64(row[j])
				x2sum[j] += uint64(row[j] * row[j])
			}
		}

		Qloc, _, _, _ = crypto.EncryptFloatMatrixRow(cryptoParams, localSketch)
		Q = mpcObj.Network.AggregateCMat(cryptoParams, Qloc)

		if debug && pid > 0 {
			pv := mpcObj.Network.CollectiveDecryptVec(cryptoParams, Q[0], 1)
			log.LLvl1("Sketch:", crypto.DecodeFloatVector(cryptoParams, pv)[:5])
			SaveMatrixToFile(cryptoParams, mpcObj, Q, nsnp, -1, pca.general.CachePath("Sketch.txt"))
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "Local bucket counts:", bucketCount)
		bucketCount = mpcObj.Network.AggregateIntVec(bucketCount)
		posCount = mpcObj.Network.AggregateIntVec(posCount)
		log.LLvl1(time.Now().Format(time.RFC3339), "Global bucket counts:", bucketCount)

		for i := range sx {
			sx[i] = rtype.FromUint64(xsum[i])
			sx2[i] = rtype.FromUint64(x2sum[i])
		}

		invN := 1.0 / float64(totInd)
		sx.MulScalar(rtype.FromFloat64(invN, 2*fracBits))
		sx2.MulScalar(rtype.FromFloat64(invN, 2*fracBits))

	}

	XMeanSS := mpcObj.TruncVec(sx, dataBits, fracBits)

	XMeanSq := mpcPar.SSSquareElemVec(XMeanSS) // E[X]^2

	sx2.Sub(XMeanSq) // E[X^2] - E[X]^2

	XVarSS := mpcObj.TruncVec(sx2, dataBits, fracBits)

	log.LLvl1(time.Now().Format(time.RFC3339), "Computing stdev ... m =", len(XVarSS))

	XStdInvSS := mpcPar.SqrtInv(XVarSS, binaryVersion)

	log.LLvl1(time.Now().Format(time.RFC3339), "Computing stdev finished")

	if pid > 0 {
		inRmat := mpc_core.InitRMat(rtype.Zero(), 3, slots*(1+((nsnp-1)/slots)))

		copy(inRmat[0], XStdInvSS)
		copy(inRmat[1], XMeanSS)
		copy(inRmat[2], XVarSS)

		outCm := mpcObj.SSToCMat(cryptoParams, inRmat)

		XStdInv = outCm[0]
		XMean = outCm[1]
		XVar = outCm[2]
	}

	if debug {
		SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{XMean}, slots*len(XMean), -1, pca.general.CachePath("XMean.txt"))
		SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{XVar}, slots*len(XVar), -1, pca.general.CachePath("XVar.txt"))
		SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{XStdInv}, slots*len(XStdInv), -1, pca.general.CachePath("XStdInv.txt"))
	}

	if !skipPowerIter {

		if restartIter <= 0 {

			if pid > 0 {
				// Normalize to reduce value range of Q
				for b := range localSketch {
					countSqrtInv := 1.0 / math.Sqrt(float64(bucketCount[b]))
					buf := make([]float64, slots)
					for i := range buf {
						buf[i] = countSqrtInv
					}
					pt, _ := crypto.EncodeFloatVector(cryptoParams, buf)

					// Cumulative weight on the mean shift (need to correct sum by meanWeight * XMean)
					meanWeight := 2*int(posCount[b]) - int(bucketCount[b])

					// Compute (Q * (1/bucketCount) - XMean) * XStdInv
					cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
						for i := range Q[b] {
							eval.MultByConstAndAdd(XMean[i], -meanWeight, Q[b][i])
							eval.MulRelin(Q[b][i], pt[0], Q[b][i])
							eval.Rescale(Q[b][i], cryptoParams.Params.Scale(), Q[b][i])
							eval.MulRelin(Q[b][i], XStdInv[i], Q[b][i])
							eval.Rescale(Q[b][i], cryptoParams.Params.Scale(), Q[b][i])
						}
						return nil
					})
				}
			} else {
				Q = make(crypto.CipherMatrix, kp)
				Qloc = make(crypto.CipherMatrix, kp)
			}

			if debug && pid > 0 {
				pv := mpcObj.Network.CollectiveDecryptVec(cryptoParams, Q[0], 1)
				log.LLvl1(time.Now().Format(time.RFC3339), "Scaling", crypto.DecodeFloatVector(cryptoParams, pv)[:5])
				SaveMatrixToFile(cryptoParams, mpcObj, Q, nsnp, -1, pca.general.CachePath("Qinit.txt"))
			}

			Q = mpcObj.Network.CollectiveBootstrapMat(cryptoParams, Q, -1)

			log.LLvl1(time.Now().Format(time.RFC3339), "Initial distributed QR, local input ", len(Q), "by", len(Q[0]), "ciphertexts")
			if pid > 0 {
				Qloc = QXLazyNormStream(cryptoParams, mpcObj, Q, XTcache, XMean, XStdInv, nRowsAll[pid])
				Qloc = mpcObj.Network.BootstrapMatAll(cryptoParams, Qloc)
				Qloc = crypto.CMultConstMat(cryptoParams, Qloc, numSnpSqrtInv, true) // scale by 1/sqrt(m)
			}

			if debug && pid > 0 {
				pv := mpcObj.Network.CollectiveDecryptVec(cryptoParams, Qloc[0], 1)
				log.LLvl1(time.Now().Format(time.RFC3339), "Before DQR", crypto.DecodeFloatVector(cryptoParams, pv)[:5])
				for outp := 1; outp < mpcObj.GetNParty(); outp++ {
					SaveMatrixToFile(cryptoParams, mpcObj, Qloc, nRowsAll[outp], outp, pca.general.CachePath("QinitX.txt"))
				}
			}

			Q = NetDQRenc(cryptoParams, mpcObj, Qloc, nRowsAll) // kp by nind

			if debug && pid > 0 {
				pv := mpcObj.Network.CollectiveDecryptVec(cryptoParams, Q[0], 1)
				log.LLvl1(time.Now().Format(time.RFC3339), "After DQR", crypto.DecodeFloatVector(cryptoParams, pv)[:5])
				for outp := 1; outp < mpcObj.GetNParty(); outp++ {
					SaveMatrixToFile(cryptoParams, mpcObj, Q, nRowsAll[outp], outp, pca.general.CachePath("QinitXOrth.txt"))
				}
			}

		} else { // Load in cached Q
			log.LLvl1(time.Now().Format(time.RFC3339), "Restarting power iteration from iter ", restartIter+1, "/", nPowerIter)

			if pid > 0 {
				// TODO cache ciphertexts instead
				cacheFile := pca.general.CachePath(fmt.Sprintf("QmulB_%d.txt", restartIter))
				mat := LoadMatrixFromFileFloat(cacheFile, ',')
				Qloc, _, _, _ = crypto.EncryptFloatMatrixRow(cryptoParams, mat)

				log.LLvl1(time.Now().Format(time.RFC3339), "Cache loaded. Number of rows:", kp, cacheFile)
			} else {
				Qloc = make(crypto.CipherMatrix, kp)
			}

			if restartIter == nPowerIter-1 {
				Q = Qloc
			} else {
				Q = NetDQRenc(cryptoParams, mpcObj, Qloc, nRowsAll)
			}

		}

		itStart := 0
		if restartIter > 0 {
			itStart = restartIter + 1
		}

		// Power iteration
		for it := itStart; it < nPowerIter; it++ {
			log.LLvl1(time.Now().Format(time.RFC3339), "Power iteration iter ", it+1, "/", nPowerIter)

			// Compute Q*X', row-based encoding
			if pid > 0 {
				Qloc := QXtLazyNormStream(cryptoParams, mpcObj, Q, Xcache, XMean, XStdInv)

				Qloc = crypto.CMultConstMat(cryptoParams, Qloc, numTotIndSqrtInv, true) // scale by 1/sqrt(n)
				Q = mpcObj.Network.AggregateCMat(cryptoParams, Qloc)
				Q = mpcObj.Network.CollectiveBootstrapMat(cryptoParams, Q, -1)
			}

			if pid > 0 {
				Qloc = QXLazyNormStream(cryptoParams, mpcObj, Q, XTcache, XMean, XStdInv, nRowsAll[pid])
				Qloc = mpcObj.Network.BootstrapMatAll(cryptoParams, Qloc)
				Qloc = crypto.CMultConstMat(cryptoParams, Qloc, numSnpSqrtInv, true) // scale by 1/sqrt(m)
			}

			if debug && pid > 0 {
				pv := mpcObj.Network.CollectiveDecryptVec(cryptoParams, Qloc[0], 1)
				log.LLvl1(time.Now().Format(time.RFC3339), "Power iter", it+1, crypto.DecodeFloatVector(cryptoParams, pv)[:5])
				for outp := 1; outp < mpcObj.GetNParty(); outp++ {
					SaveMatrixToFile(cryptoParams, mpcObj, Qloc, nRowsAll[outp], outp, pca.general.CachePath(fmt.Sprintf("QmulB_%d.txt", it)))
				}
			}

			// Skip QR in the last iteration
			if it == nPowerIter-1 {
				Q = Qloc
			} else {
				Q = NetDQRenc(cryptoParams, mpcObj, Qloc, nRowsAll)
			}
		}
		log.LLvl1(time.Now().Format(time.RFC3339), "Power iteration complete")

		if debug && pid > 0 {
			pv := mpcObj.Network.CollectiveDecryptVec(cryptoParams, Q[0], 1)
			log.LLvl1(time.Now().Format(time.RFC3339), "After power iter", crypto.DecodeFloatVector(cryptoParams, pv)[:5])
			for outp := 1; outp < mpcObj.GetNParty(); outp++ {
				SaveMatrixToFile(cryptoParams, mpcObj, Q, nRowsAll[outp], outp, pca.general.CachePath(fmt.Sprintf("Q_final.txt")))
			}
		}

	} else {
		log.LLvl1(time.Now().Format(time.RFC3339), "Power iteration skipped. Using Q_final from a previous run.")

		if pid > 0 {
			// TODO cache ciphertexts instead
			mat := LoadMatrixFromFileFloat(pca.general.CachePath("QmulB_9.txt"), ',') // TODO: Temporary fix, fetch Q_final instead
			Q, _, _, _ = crypto.EncryptFloatMatrixRow(cryptoParams, mat)
		} else {
			Q = make(crypto.CipherMatrix, kp)
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "Cache loaded. Number of rows:", kp)
	}

	// Q contains Q*X' (kp by numInd) for each party
	// Compute local Gram matrix Q * Q'
	// TODO: be careful of the increasing data range
	if pid > 0 {

		log.LLvl1(time.Now().Format(time.RFC3339), "Computing covariance matrix")

		nct := ((kp*kp)-1)/slots + 1
		Zloc := crypto.CZeros(cryptoParams, nct)
		for i := 0; i < kp; i++ {
			for j := i; j < kp; j++ {
				inds := [2]int{i*kp + j, j*kp + i}

				iprod := crypto.InnerProd(cryptoParams, Q[i], Q[j])

				for k := range inds {
					ctid, slotid := inds[k]/slots, inds[k]%slots

					ct := crypto.Mask(cryptoParams, iprod, slotid, false)

					cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
						eval.Add(ct, Zloc[ctid], Zloc[ctid])
						return nil
					})

					if i == j {
						break
					}
				}
			}
		}

		Z := mpcObj.Network.AggregateCVec(cryptoParams, Zloc)

		// Normalize by 1/N
		Z = crypto.CMultConst(cryptoParams, Z, 1.0/float64(totInd), true)

		if debug {
			SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{Z}, kp*kp, -1, pca.general.CachePath("Zgram.txt"))
		}

		Zss := mpcObj.CVecToSS(cryptoParams, mpcObj.GetRType(), Z, -1, len(Z), kp*kp)

		for i := range Zmat {
			Zmat[i] = Zss[(i * kp):((i + 1) * kp)]
		}
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Eigen decomposition")

	// Eigen decomposition
	Vss, L := mpcObj.EigenDecomp(Zmat)
	Vss, L = mpcObj.SortRowsDescend(Vss, L)
	Vss = Vss[:npc]

	if debug && pid > 0 {
		Vr := mpcObj.RevealSymMat(Vss)
		Lr := mpcObj.RevealSymVec(L)
		Vf := Vr.ToFloat(fracBits)
		for i := range Vf {
			log.LLvl1(time.Now().Format(time.RFC3339), "V[i]", i, Vf[i])
		}
		log.LLvl1(time.Now().Format(time.RFC3339), "L", Lr.ToFloat(fracBits))
	}

	V := mpcObj.SSToCMat(cryptoParams, Vss)

	if debug && pid > 0 {
		SaveMatrixToFile(cryptoParams, mpcObj, V, kp, -1, pca.general.CachePath("V.txt"))
	}

	Qpc := crypto.CZeroMat(cryptoParams, len(Q[0]), npc)
	if pid > 0 {
		log.LLvl1(time.Now().Format(time.RFC3339), "Extract PC subspace")

		// Extract local PC subspace by computing V*Q (npc by numInd)
		for r := range V {
			for c := range Q {
				ctid, slotid := c/slots, c%slots

				elem := crypto.Mask(cryptoParams, V[r][ctid], slotid, false)
				elem = crypto.InnerSumAll(cryptoParams, crypto.CipherVector{elem})

				cv := crypto.CMultScalar(cryptoParams, Q[c], elem)

				cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
					for i := range cv {
						eval.Add(cv[i], Qpc[r][i], Qpc[r][i])
					}
					return nil
				})
			}
		}
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "AssertSync")
	mpcObj.AssertSync()

	return Qpc
}
