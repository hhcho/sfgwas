package gwas

import (
	"fmt"
	"math"
	"runtime"
	"sync"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/ldsec/lattigo/v2/ckks"
	"go.dedis.ch/onet/v3/log"

	"github.com/hhcho/sfgwas/crypto"

	"gonum.org/v1/gonum/mat"
)

type AssocTest struct {
	general *ProtocolInfo

	pheno crypto.PlainVector

	inputCov crypto.PlainMatrix
	Qpc      crypto.CipherMatrix

	MafLowerBound float64
	MafUpperBound float64
}

func (g *ProtocolInfo) InitAssociationTests(Qpc crypto.CipherMatrix) *AssocTest {

	pid := g.mpcObj[0].GetPid()
	gwasParams := g.gwasParams
	cps := g.cps

	var phenoEnc crypto.PlainVector
	var covEnc crypto.PlainMatrix

	if pid > 0 {
		phenoEnc = crypto.EncodeDense(cps, mat.DenseCopyOf(g.pheno))[0]
		covEnc = crypto.EncodeDense(cps, mat.DenseCopyOf(g.cov))
		r, c := g.cov.Dims()
		log.LLvl1(time.Now().Format(time.RFC3339), "Cov dims:", r, c)
	} else {
		phenoEnc = make(crypto.PlainVector, 0)
		covEnc = make(crypto.PlainMatrix, gwasParams.NumCov())
		log.LLvl1(time.Now().Format(time.RFC3339), "Cov dims:", 0, gwasParams.NumCov())
	}

	return &AssocTest{
		general:  g,
		pheno:    phenoEnc,
		inputCov: covEnc,
		Qpc:      Qpc,
	}
}

// To add for UKB
// func (ast *AssocTest) Preprocessing() {
// }

// func (ast *AssocTest) AssociationTest(blockIndex, start, end int) (crypto.CipherVector, []bool) {
// }

// Orthogonal basis of covariates and PCs combined; joint QR
func (ast *AssocTest) computeCombinedQV2(C crypto.PlainMatrix, Qpc crypto.CipherMatrix) crypto.CipherMatrix {
	cryptoParams := ast.general.cps
	mpcPar := ast.general.mpcObj
	mpcObj := mpcPar[0]
	pid := mpcPar[0].GetPid()
	slots := cryptoParams.GetSlots()

	gwasParams := ast.general.gwasParams
	nrowsAll := gwasParams.FiltNumInds()
	nrowsTotal := 0
	for i := 1; i < len(nrowsAll); i++ {
		nrowsTotal += nrowsAll[i]
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Starting QR: C numCols", len(C))

	CEnc := crypto.EncryptPlaintextMatrix(cryptoParams, C)

	comb := make(crypto.CipherMatrix, len(C)+len(Qpc))
	for i := range CEnc {
		comb[i] = CEnc[i]
	}
	for i := range Qpc {
		comb[len(CEnc)+i] = Qpc[i]
	}

	start := time.Now()

	Qcomb := NetDQRenc(cryptoParams, mpcObj, comb, nrowsAll)

	log.LLvl1(time.Now().Format(time.RFC3339), "Covariate joint QR time: ", time.Since(start))
	log.LLvl1(time.Now().Format(time.RFC3339), "Qcomb dimensions: r,c :", len(Qcomb[0]), len(Qcomb))

	Qcomb = mpcObj.Network.BootstrapMatAll(cryptoParams, Qcomb)

	log.LLvl1(time.Now().Format(time.RFC3339), "Qcomb replacing first vector with an all-ones vector (normalized)")
	if pid > 0 {
		ct := crypto.CZeros(cryptoParams, 1)[0]
		ct = crypto.AddConst(cryptoParams, ct, 1.0)

		QFirst := make(crypto.CipherVector, ((nrowsAll[pid]-1)/slots)+1)
		for i := range QFirst {
			nElem := slots
			if i == len(QFirst)-1 {
				nElem = nrowsAll[pid] - (len(QFirst)-1)*slots
			}
			QFirst[i] = crypto.MaskTrunc(cryptoParams, ct, nElem)
		}

		Qcomb[0] = QFirst
		Qcomb, _ = crypto.FlattenLevels(cryptoParams, Qcomb)
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "AssertSync")
	mpcObj.AssertSync()

	return Qcomb
}

// Orthogonal basis of covariates and PCs combined
func (ast *AssocTest) computeCombinedQ(C crypto.PlainMatrix, Qpc crypto.CipherMatrix) crypto.CipherMatrix {
	debug := ast.general.config.Debug
	cryptoParams := ast.general.cps
	mpcPar := ast.general.mpcObj
	mpcObj := mpcPar[0]
	pid := mpcPar[0].GetPid()
	slots := cryptoParams.GetSlots()

	gwasParams := ast.general.gwasParams
	nrowsAll := gwasParams.FiltNumInds()
	nrowsTotal := 0
	for i := 1; i < len(nrowsAll); i++ {
		nrowsTotal += nrowsAll[i]
	}
	nrowsTotalInv := 1.0 / float64(nrowsTotal)

	ncov := gwasParams.NumCov()
	npc := gwasParams.NumPC()

	log.LLvl1(time.Now().Format(time.RFC3339), "Starting QR: C numCols", len(C))

	start := time.Now()

	// QR on covariates
	var Qi crypto.CipherMatrix

	CEnc := crypto.EncryptPlaintextMatrix(cryptoParams, C)
	Qi = NetDQRenc(cryptoParams, mpcObj, CEnc, nrowsAll)

	// Alternative approach (DASH): appears to be less accurate, needs a closer look
	// 	Qi = NetDQRplain(cryptoParams, mpcObj, C, nrowsAll)

	log.LLvl1(time.Now().Format(time.RFC3339), "Covariate QR time: ", time.Since(start))
	log.LLvl1(time.Now().Format(time.RFC3339), "Qi dimensions: r,c :", len(Qi[0]), len(Qi))

	Qi = mpcObj.Network.BootstrapMatAll(cryptoParams, Qi)

	log.LLvl1(time.Now().Format(time.RFC3339), "Qi replacing first vector with an all-ones vector (normalized)")
	if pid > 0 {
		ct := crypto.CZeros(cryptoParams, 1)[0]
		ct = crypto.AddConst(cryptoParams, ct, 1.0)

		QiFirst := make(crypto.CipherVector, ((nrowsAll[pid]-1)/slots)+1)
		for i := range QiFirst {
			nElem := slots
			if i == len(QiFirst)-1 {
				nElem = nrowsAll[pid] - (len(QiFirst)-1)*slots
			}
			QiFirst[i] = crypto.MaskTrunc(cryptoParams, ct, nElem)
		}

		Qi[0] = QiFirst
		Qi, _ = crypto.FlattenLevels(cryptoParams, Qi)
	}

	if debug && pid > 0 {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			SaveMatrixToFile(cryptoParams, mpcObj, Qi, nrowsAll[party], party, ast.general.CachePath("Qi.txt"))
		}
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "AssertSync")
	mpcObj.AssertSync()

	// Compute (I - Qi * Qi') * Qpc
	fn := func(cp *crypto.CryptoParams, a crypto.CipherVector,
		B crypto.CipherMatrix, j int) crypto.CipherVector {
		return crypto.CMult(cp, a, B[j])
	}

	Qpcmi := make(crypto.CipherMatrix, npc)

	if pid > 0 && Qpc != nil {

		Qpcmi = DCMatMulAAtB(cryptoParams, mpcObj, Qi, Qpc, nrowsAll, npc, fn)

		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				SaveMatrixToFile(cryptoParams, mpcObj, Qpc, nrowsAll[party], party, ast.general.CachePath("Qpc_assoc.txt"))
				SaveMatrixToFile(cryptoParams, mpcObj, Qpcmi, nrowsAll[party], party, ast.general.CachePath("Qpcmi_mult.txt"))
			}
		}

		for c := range Qpcmi {
			Qpcmi[c] = crypto.CMultConstRescale(cryptoParams, Qpcmi[c], nrowsTotalInv, true)
			Qpcmi[c] = crypto.CSub(cryptoParams, Qpc[c], Qpcmi[c])
		}

		Qpcmi = mpcObj.Network.BootstrapMatAll(cryptoParams, Qpcmi)
	}

	if debug && pid > 0 {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			SaveMatrixToFile(cryptoParams, mpcObj, Qpcmi, nrowsAll[party], party, ast.general.CachePath("Qpcmi_before_QR.txt"))
		}
	}

	// QR factorize Qpcmi
	if Qpc != nil {
		start = time.Now()

		Qpcmi = NetDQRenc(cryptoParams, mpcObj, Qpcmi, nrowsAll)

		log.LLvl1(time.Now().Format(time.RFC3339), "QR(Qpcmi) time: ", time.Since(start))
		log.LLvl1(time.Now().Format(time.RFC3339), "QR(Qpcmi) dimensions: r,c :", len(Qi[0]), len(Qi))

		Qpcmi = mpcObj.Network.BootstrapMatAll(cryptoParams, Qpcmi)
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "AssertSync")
	mpcObj.AssertSync()

	if debug && pid > 0 {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			SaveMatrixToFile(cryptoParams, mpcObj, Qpcmi, nrowsAll[party], party, ast.general.CachePath("Qpcmi.txt"))
		}
	}

	var Q crypto.CipherMatrix

	if pid > 0 {
		// Concatenate Qi and Qpcmi
		if Qpc == nil {
			Q = Qi
		} else {
			Q = make(crypto.CipherMatrix, ncov+npc)
			for c := range Qi {
				Q[c] = Qi[c]
			}
			for c := range Qpcmi {
				Q[ncov+c] = Qpcmi[c]
			}
		}
	}

	return Q
}

func (ast *AssocTest) GenoBlockMult(b int, mat crypto.CipherMatrix) (matOut crypto.CipherMatrix, dosageSum, dosageSqSum []float64, filtOut []bool) {
	cryptoParams := ast.general.cps

	pid := ast.general.mpcObj[0].GetPid()

	slots := cryptoParams.GetSlots()
	isPgen := ast.general.IsPgen()
	pgenBatchSize := ast.general.config.PgenBatchSize

	XBlock := ast.general.genoBlocks[b]
	numBlocks := ast.general.config.GenoNumBlocks

	gwasParams := ast.general.gwasParams
	snpFilt := gwasParams.snpFilt

	blockSize := ast.general.genoBlockSizes[b]

	shift := uint64(0)
	for i := 0; i < b; i++ {
		shift += uint64(ast.general.genoBlockSizes[i])
	}

	var nsnps int
	if isPgen {
		if snpFilt == nil { // no QC
			nsnps = blockSize
		} else {
			nsnps = SumBool(snpFilt[shift : shift+uint64(blockSize)])
		}
	} else {
		nsnps = int(XBlock.NumColsToKeep())
	}

	if nsnps == 0 { // empty block
		log.LLvl1(time.Now().Format(time.RFC3339), "MatMult: block", b+1, "/", numBlocks, "skipped (empty)")
		return
	}

	numCtx := 0
	if isPgen {
		for nleft := nsnps; nleft > 0; {
			bsize := Min(nleft, pgenBatchSize)
			numCtx += 1 + (bsize-1)/slots
			nleft -= bsize
		}
	} else {
		numCtx = 1 + (nsnps-1)/slots
	}

	multFile := ast.general.CachePath(fmt.Sprintf("assoc_cache_mult.%d.bin", b))
	dosFile := ast.general.CachePath(fmt.Sprintf("assoc_cache_dos_sum.%d.txt", b))
	dos2File := ast.general.CachePath(fmt.Sprintf("assoc_cache_dos_sqsum.%d.txt", b))
	filtFile := ast.general.CachePath(fmt.Sprintf("assoc_cache_filt.%d.txt", b))

	if fileExists(multFile) && fileExists(dosFile) && fileExists(dos2File) && fileExists(filtFile) {

		log.LLvl1(time.Now().Format(time.RFC3339), "MatMult: block", b+1, "/", numBlocks, "cache found")

		matOut = crypto.LoadCipherMatrixFromFile(cryptoParams, multFile)
		dosageSum = LoadFloatVectorFromFile(dosFile, numCtx*slots)
		dosageSqSum = LoadFloatVectorFromFile(dos2File, numCtx*slots)
		filtOut = readFilterFromFile(filtFile, numCtx*slots, true)

		log.LLvl1("Dosage Sum:", dosageSum[:5])
		log.LLvl1("Dosage SqSum:", dosageSqSum[:5])

	} else {

		log.LLvl1(time.Now().Format(time.RFC3339), "MatMult: block", b+1, "/", numBlocks, "starting")

		filtOut = make([]bool, numCtx*slots)

		start := time.Now()
		if isPgen {

			numInd := gwasParams.numFiltInds[pid]
			snpFilt := gwasParams.snpFilt[shift : shift+uint64(blockSize)]

			pgenFile := fmt.Sprintf(ast.general.config.GenoFilePrefix, b+1) // Geno file for chromosome b+1

			startIndex := 0
			counter := 0
			outShift := 0
			batchIndex := 0

			nbatch := 1 + (nsnps-1)/pgenBatchSize
			outMult := make([]crypto.CipherMatrix, nbatch)

			dosageSum = make([]float64, numCtx*slots)
			dosageSqSum = make([]float64, numCtx*slots)

			nblocksParallel := ast.general.config.LocalAssocNumBlocksParallel
			nprocsPerBlock := Max(1, runtime.GOMAXPROCS(0)/nblocksParallel)

			// Dispatcher
			threadPool := make(chan int, nblocksParallel)
			for i := 0; i < nblocksParallel; i++ {
				threadPool <- i
			}

			var wg sync.WaitGroup

			for idx := 0; idx < blockSize; idx++ {
				if snpFilt[idx] {
					counter++
				}

				if counter == pgenBatchSize || (idx == blockSize-1 && counter > 0) {
					// Fetch an available thread
					threadId := <-threadPool
					wg.Add(1)

					go func(threadId, batchIndex, startIndex, idx, counter, shift, outShift int) {
						defer wg.Done()

						start := time.Now()
						log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("MatMult: block %d/%d, batch %d/%d, thread %d started", b+1, numBlocks, batchIndex+1, nbatch, threadId))

						batchFilt := snpFilt[startIndex : idx+1]
						gfsTempFile := ast.general.CachePath(fmt.Sprintf("pgen_gfs.%d.tmp", threadId))
						FilterMatrixFilePgen(pgenFile, numInd, counter, ast.general.config.SampleKeepFile, ast.general.config.SnpIdsFile, shift+startIndex, batchFilt, gfsTempFile)

						X := NewGenoFileStream(gfsTempFile, uint64(numInd), uint64(counter), true)

						mult, sum, sqSum := MatMult4Stream(cryptoParams, mat, X, 5, true, nprocsPerBlock)

						outMult[batchIndex] = mult
						copy(dosageSum[outShift:outShift+len(sum)], sum)
						copy(dosageSqSum[outShift:outShift+len(sqSum)], sqSum)
						for c := 0; c < counter; c++ {
							filtOut[outShift+c] = true
						}

						log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("MatMult: block %d/%d, batch %d/%d, thread %d finished,", b+1, numBlocks, batchIndex, nbatch, threadId), "elapsed time", time.Since(start))

						// Return thread to pool
						threadPool <- threadId
					}(threadId, batchIndex, startIndex, idx, counter, int(shift), outShift)

					nctx := 1 + (counter-1)/slots
					outShift += nctx * slots
					startIndex = idx + 1
					batchIndex++
					counter = 0
				}
			}

			wg.Wait() // Wait until all batches are finished
			close(threadPool)

			matOut = crypto.ConcatCipherMatrix(outMult)

		} else {
			matOut, dosageSum, dosageSqSum = MatMult4Stream(cryptoParams, mat, XBlock, 5, true, 0)

			for c := 0; c < nsnps; c++ {
				filtOut[c] = true
			}
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "MatMult: block", b+1, "/", numBlocks, "elapsed time", time.Since(start))

		// Save cache
		crypto.SaveCipherMatrixToFile(cryptoParams, matOut, multFile)
		SaveFloatVectorToFile(dosFile, dosageSum)
		SaveFloatVectorToFile(dos2File, dosageSqSum)
		writeFilterToFile(filtFile, filtOut, true)
	}

	return
}

func (ast *AssocTest) GetAssociationStats() (crypto.CipherVector, []bool) {
	debug := ast.general.config.Debug

	covAllOnes := ast.general.config.CovAllOnes // Flag indicating whether cov includes an all-ones covariate

	cryptoParams := ast.general.cps
	mpcPar := ast.general.mpcObj
	mpcObj := mpcPar[0]
	pid := mpcPar[0].GetPid()
	gwasParams := ast.general.gwasParams
	slots := cryptoParams.GetSlots()

	numBlocks := ast.general.config.GenoNumBlocks

	/* Sample counts */
	nrowsAll := gwasParams.FiltNumInds()
	nrowsTotal := 0
	for i := 1; i < len(nrowsAll); i++ {
		nrowsTotal += nrowsAll[i]
	}
	nrowsTotalInv := 1.0 / float64(nrowsTotal)

	/* Phenotypes and PCs */
	y := ast.pheno
	Qpc := ast.Qpc

	/* Setup covariates */
	C := ast.inputCov
	ncov := gwasParams.NumCov()
	if !covAllOnes {
		log.LLvl1("Adding an all-ones covariate")

		arr := make([]float64, nrowsAll[pid])
		for i := range arr {
			arr[i] = 1.0
		}
		pv, _ := crypto.EncodeFloatVector(cryptoParams, arr)

		C = append([]crypto.PlainVector{pv}, C...)
		ncov += 1

		covAllOnes = true
	} else {
		log.LLvl1("Warning: assumes the first covariate is all ones (if not reorder)")
	}

	if debug && pid > 0 {
		yf := crypto.DecodeFloatVector(cryptoParams, y)[:nrowsAll[pid]]
		SaveFloatVectorToFile(ast.general.CachePath("y.txt"), yf)

		Cf := make([][]float64, len(C))
		for i := range C {
			Cf[i] = crypto.DecodeFloatVector(cryptoParams, C[i])[:nrowsAll[pid]]
		}
		SaveFloatMatrixToFile(ast.general.CachePath("C.txt"), Cf)
	}

	cacheFileQ := ast.general.CachePath("Qcomb.bin")
	var Q crypto.CipherMatrix
	if ast.general.config.UseCachedCombinedQ {
		if pid > 0 {
			Q = crypto.LoadCipherMatrixFromFile(cryptoParams, cacheFileQ)
			log.LLvl1(time.Now().Format(time.RFC3339), "Qcomb loaded from", cacheFileQ)
		}
	} else {
		Q = ast.computeCombinedQV2(C, Qpc) // nil for pid = 0
		if pid > 0 {
			crypto.SaveCipherMatrixToFile(cryptoParams, Q, cacheFileQ)
			log.LLvl1(time.Now().Format(time.RFC3339), "Qcomb saved to", cacheFileQ)
		}
	}

	if debug && pid > 0 {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			SaveMatrixToFile(cryptoParams, mpcObj, Q, nrowsAll[party], party, ast.general.CachePath("Qcomb.txt"))
		}
	}

	var varx, vary, sx, sy, sxy crypto.CipherVector
	var nsnps, numCtx int
	var outFilter []bool

	if pid == 0 {
		numCtx = mpcObj.Network.ReceiveInt(mpcObj.GetHubPid())
		nsnps = mpcObj.Network.ReceiveInt(mpcObj.GetHubPid())

		varx = crypto.CZeros(cryptoParams, numCtx)
		vary = crypto.CZeros(cryptoParams, 1)

	} else { // pid > 0

		// Project covariates out of y: ynew = (I - Q*Q')*y
		ymat := make(crypto.PlainMatrix, 1)
		ymat[0] = y

		mmplainfn := func(cp *crypto.CryptoParams, a crypto.CipherVector,
			B crypto.PlainMatrix, j int) crypto.CipherVector {
			return crypto.CPMult(cp, a, B[j])
		}

		ynew := DCMatMulAAtBPlain(cryptoParams, mpcObj, Q, ymat, nrowsAll, 1, mmplainfn) // Level -2
		ynew[0] = crypto.CMultConstRescale(cryptoParams, ynew[0], nrowsTotalInv, true)

		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				SaveMatrixToFile(cryptoParams, mpcObj, ynew, nrowsAll[party], party, ast.general.CachePath("QQy.txt"))
			}
		}

		ynew[0] = mpcObj.Network.BootstrapVecAll(cryptoParams, ynew[0])
		ynew[0] = crypto.CMultConst(cryptoParams, ynew[0], -1.0, true)
		ynew[0] = crypto.CPAdd(cryptoParams, ynew[0], y)

		log.LLvl1(time.Now().Format(time.RFC3339), "ynew computed")

		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				SaveMatrixToFile(cryptoParams, mpcObj, ynew, nrowsAll[party], party, ast.general.CachePath("ynew.txt"))
			}
		}

		// Compute u (=Q*Q'*1)
		dummyMat := make(crypto.PlainMatrix, 1)
		mm1fn := func(cp *crypto.CryptoParams, a crypto.CipherVector,
			B crypto.PlainMatrix, j int) crypto.CipherVector {
			return crypto.CopyEncryptedVector(a)
		}

		u := DCMatMulAAtBPlain(cryptoParams, mpcObj, Q, dummyMat, nrowsAll, 1, mm1fn) // Level -2
		u[0] = crypto.CMultConstRescale(cryptoParams, u[0], nrowsTotalInv, true)
		log.LLvl1(time.Now().Format(time.RFC3339), "u computed")

		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				SaveMatrixToFile(cryptoParams, mpcObj, u, nrowsAll[party], party, ast.general.CachePath("u.txt"))
			}
		}

		omu := crypto.CZeros(cryptoParams, len(u[0]))
		if !covAllOnes {
			omu = crypto.CSub(cryptoParams, omu, u[0])
			omu = crypto.CAddConst(cryptoParams, omu, 1.0)

			log.LLvl1(time.Now().Format(time.RFC3339), "omu computed")
		} else {
			log.LLvl1(time.Now().Format(time.RFC3339), "omu set to zero")
		}

		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{omu}, nrowsAll[party], party, ast.general.CachePath("omu.txt"))
			}
		}

		// Compute sx, sxx, sxy, sy, syy
		// sx = 1'*X - u'*X = omu'*X
		// sxx = diag(X'X) - diag(B'B)
		// sxy = ynew'*X

		// In parallel:
		// (1) Compute B (=Q'*X, row-based encoding)
		// (2) Compute sx (=omu'*X)
		// (3) Compute sxy (=ynew'*X)
		// Note: if covAllOnes = true, then sx = sy = 0. Skip all calculations involving sx and sy.

		concat := make(crypto.CipherMatrix, len(Q)+2) // remove all ones
		for i := 0; i < len(Q); i++ {
			concat[i] = Q[i]
		}
		concat[len(Q)] = omu
		concat[len(Q)+1] = ynew[0]

		filtOut := make([][]bool, numBlocks)

		log.LLvl1(time.Now().Format(time.RFC3339), "Multiplication with genotype matrix started")

		sxBlocks := make([]crypto.CipherMatrix, numBlocks)
		sxxBlocks := make([]crypto.CipherMatrix, numBlocks)
		sxyBlocks := make([]crypto.CipherMatrix, numBlocks)

		for b := 0; b < numBlocks; b++ {
			if !ast.general.IsBlockForAssocTest(b) {
				log.LLvl1(time.Now().Format(time.RFC3339), "MatMult: block", b+1, "/", numBlocks, "skipped")
			} else {
				concatOut, dosageSum, dosageSqSum, filt := ast.GenoBlockMult(b, concat)
				if concatOut == nil {
					continue
				}

				log.LLvl1(time.Now().Format(time.RFC3339), "block", b+1, "/", numBlocks, "aggregating")
				concatOut = mpcObj.Network.AggregateCMat(cryptoParams, concatOut)
				log.LLvl1(time.Now().Format(time.RFC3339), "block", b+1, "/", numBlocks, "bootstrapping")
				concatOut = mpcObj.Network.CollectiveBootstrapMat(cryptoParams, concatOut, -1)

				B := make(crypto.CipherMatrix, len(Q)-1) // Skip the one correponding to all ones
				for i := range B {
					B[i] = crypto.CMultConstRescale(cryptoParams, concatOut[i+1], math.Sqrt(nrowsTotalInv), true)
				}

				if covAllOnes {
					sxBlocks[b] = crypto.CipherMatrix{crypto.CZeros(cryptoParams, len(concatOut[len(Q)]))}
					log.LLvl1(time.Now().Format(time.RFC3339), "sx set to zero")
				} else {
					sxBlocks[b] = crypto.CipherMatrix{concatOut[len(Q)]}
				}

				sxyBlocks[b] = crypto.CipherMatrix{concatOut[len(Q)+1]}

				log.LLvl1(time.Now().Format(time.RFC3339), "block", b+1, "/", numBlocks, "computed B, sx, sxy")

				var sx2 crypto.CipherVector
				if dosageSqSum != nil {
					sxxBlocks[b] = make(crypto.CipherMatrix, 1)
					sxxBlocks[b][0], _ = crypto.EncryptFloatVector(cryptoParams, dosageSqSum)
					sx2, _ = crypto.EncryptFloatVector(cryptoParams, dosageSum)
				}

				sx2 = mpcObj.Network.AggregateCVec(cryptoParams, sx2)
				sx2 = crypto.CMultConstRescale(cryptoParams, sx2, math.Sqrt(nrowsTotalInv), true)

				if pid == mpcObj.GetHubPid() {
					cryptoParams.WithEvaluator(func(evaluator ckks.Evaluator) error {
						for c := range B {
							for j := range sxxBlocks[b][0] {
								tmp := evaluator.MulRelinNew(B[c][j], B[c][j])
								evaluator.Sub(sxxBlocks[b][0][j], tmp, sxxBlocks[b][0][j])
							}
						}
						for j := range sxxBlocks[b][0] {
							tmp := evaluator.MulRelinNew(sx2[j], sx2[j])
							evaluator.Sub(sxxBlocks[b][0][j], tmp, sxxBlocks[b][0][j])
						}
						return nil
					})
				}

				sxxBlocks[b][0] = mpcObj.Network.AggregateCVec(cryptoParams, sxxBlocks[b][0])

				log.LLvl1(time.Now().Format(time.RFC3339), "block", b+1, "/", numBlocks, "computed sxx")

				filtOut[b] = filt
			}
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "All blocks processed")

		sx = crypto.ConcatCipherMatrix(sxBlocks)[0]
		sxy = crypto.ConcatCipherMatrix(sxyBlocks)[0]
		sxx := crypto.ConcatCipherMatrix(sxxBlocks)[0]
		sxBlocks, sxxBlocks, sxyBlocks = nil, nil, nil

		totLen := 0
		for i := range filtOut {
			totLen += len(filtOut[i])
		}

		outFilter = make([]bool, totLen)
		outShift := 0
		for i := range filtOut {
			copy(outFilter[outShift:], filtOut[i])
			outShift += len(filtOut[i])
		}
		filtOut = nil

		numCtx = len(sx)
		nsnps = SumBool(outFilter)

		if pid == mpcObj.GetHubPid() {
			mpcObj.Network.SendInt(numCtx, 0)
			mpcObj.Network.SendInt(nsnps, 0)
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "numCtx", numCtx, "numSnps", nsnps)

		// Compute sy and syy
		if covAllOnes {
			sy = crypto.CZeros(cryptoParams, 1)
			log.LLvl1(time.Now().Format(time.RFC3339), "sy set to zero")
		} else {
			syloc := crypto.InnerSumAll(cryptoParams, ynew[0])
			sy = crypto.CipherVector{mpcObj.Network.AggregateCText(cryptoParams, syloc)}
			sy = mpcObj.Network.CollectiveBootstrapVec(cryptoParams, sy, -1)
		}

		ynewsq := crypto.CMult(cryptoParams, ynew[0], ynew[0])
		syyloc := crypto.InnerSumAll(cryptoParams, ynewsq)

		syy := crypto.CipherVector{mpcObj.Network.AggregateCText(cryptoParams, syyloc)}
		syy = mpcObj.Network.CollectiveBootstrapVec(cryptoParams, syy, -1)

		log.LLvl1(time.Now().Format(time.RFC3339), "Computed sy/syy")

		totalInds := 0
		for _, v := range nrowsAll {
			totalInds += v
		}

		sqrtinvn := 1.0 / math.Sqrt(float64(totalInds))
		if !covAllOnes {
			sx = crypto.CMultConst(cryptoParams, sx, sqrtinvn, true) // sx / sqrt(n)
			sy = crypto.CMultConst(cryptoParams, sy, sqrtinvn, true) // sy / sqrt(n)

			varx = crypto.CMult(cryptoParams, sx, sx)   // sx * sx / n
			varx = crypto.CSub(cryptoParams, sxx, varx) // varx = sxx - (sx * sx / n)

			vary = crypto.CMult(cryptoParams, sy, sy)   // sy * sy / n
			vary = crypto.CSub(cryptoParams, syy, vary) // vary = syy - (sy * sy / n)
		} else {
			varx = sxx
			vary = syy
		}

		if debug {
			writeFilterToFile(ast.general.CachePath("xfilt.bin"), outFilter, true)
			SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{sx}, len(sx)*slots, -1, ast.general.CachePath("sx.txt"))     // sx / sqrt(n)
			SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{sxx}, len(sx)*slots, -1, ast.general.CachePath("sxx.txt"))   // sxx
			SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{sy}, 1, -1, ast.general.CachePath("sy.txt"))                 // sy / sqrt(n)
			SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{syy}, 1, -1, ast.general.CachePath("syy.txt"))               // syy
			SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{sxy}, len(sx)*slots, -1, ast.general.CachePath("sxy.txt"))   // sxy
			SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{varx}, len(sx)*slots, -1, ast.general.CachePath("varx.txt")) // sxx - (sx*sx/n)
			SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{vary}, 1, -1, ast.general.CachePath("vary.txt"))             // syy - (sy*sy/n)
		}
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "AssertSync")
	mpcObj.AssertSync()

	stdinvx, stdinvy := ast.computeStdInv(varx, vary, nsnps, outFilter)
	log.LLvl1(time.Now().Format(time.RFC3339), "Computed stdev")

	if pid > 0 {
		var stats crypto.CipherVector
		if !covAllOnes {
			stats = crypto.CMultScalar(cryptoParams, sx, sy[0]) // sx * sy / n
			stats = crypto.CSub(cryptoParams, sxy, stats)       // sxy - (sx * sy / n)
		} else {
			stats = sxy
		}
		stats = crypto.CMult(cryptoParams, stats, stdinvx)       // stdinvx * (sxy - (sx * sy) / n)
		stats = crypto.CMultScalar(cryptoParams, stats, stdinvy) // stdinvx * stdinvy * (sxy - (sx * sy) / n)

		log.LLvl1(time.Now().Format(time.RFC3339), "All done!")

		return stats, outFilter
	}

	return nil, nil // party 0
}

// Returns stdinvx and stdinvy
func (ast *AssocTest) computeStdInv(varx, vary crypto.CipherVector, nsnps int, filter []bool) (crypto.CipherVector, *ckks.Ciphertext) {
	debug := ast.general.config.Debug

	cryptoParams := ast.general.cps
	mpcPar := ast.general.mpcObj
	mpcObj := mpcPar[0]
	pid := mpcPar[0].GetPid()
	rtype := mpcPar[0].GetRType()
	slots := cryptoParams.GetSlots()
	useBoolean := mpcPar[0].GetBooleanShareFlag()

	// Convert to SS
	varxSS := mpcObj.CVecToSS(cryptoParams, mpcObj.GetRType(), varx, -1, len(varx), slots*len(varx))
	varySS := mpcObj.CiphertextToSS(cryptoParams, mpcObj.GetRType(), vary[0], -1, 1)

	if debug && pid > 0 {
		log.LLvl1(time.Now().Format(time.RFC3339), "varxSS", mpcObj.RevealSymVec(varxSS[:5]).ToFloat(mpcObj.GetFracBits()))
		log.LLvl1(time.Now().Format(time.RFC3339), "varySS", mpcObj.RevealSymVec(varySS).ToFloat(mpcObj.GetFracBits()))
	}

	// Concatenate
	varSS := mpc_core.InitRVec(rtype.Zero(), nsnps+1)
	if pid > 0 {
		dst := 0
		for src := range varxSS {
			if filter[src] {
				varSS[dst] = varxSS[src]
				dst++
			}
		}
	}

	varSS[len(varSS)-1] = varySS[0]

	// Compute Sqrt Inverse
	stdinvSS := mpcPar.SqrtInv(varSS, useBoolean)

	if debug && pid > 0 {
		log.LLvl1(time.Now().Format(time.RFC3339), "varxSS", mpcObj.RevealSymVec(varxSS[:5]).ToFloat(mpcObj.GetFracBits()))
		log.LLvl1(time.Now().Format(time.RFC3339), "varSS", mpcObj.RevealSymVec(varSS[:5]).ToFloat(mpcObj.GetFracBits()))
		log.LLvl1(time.Now().Format(time.RFC3339), "stdinvxSS", mpcObj.RevealSymVec(stdinvSS[:5]).ToFloat(mpcObj.GetFracBits()))
		log.LLvl1(time.Now().Format(time.RFC3339), "stdinvySS", mpcObj.RevealSymVec(stdinvSS[(len(stdinvSS)-1):]).ToFloat(mpcObj.GetFracBits()))
	}

	// Convert back to HE
	stdinvxSS := mpc_core.InitRVec(rtype.Zero(), len(varxSS))
	if pid > 0 {
		src := 0
		for dst := range filter {
			if filter[dst] {
				stdinvxSS[dst] = stdinvSS[src]
				src++
			}
		}
	}

	stdinvx := mpcObj.SSToCVec(cryptoParams, stdinvxSS)
	stdinvy := mpcObj.SStoCiphertext(cryptoParams, mpc_core.RVec{stdinvSS[len(stdinvSS)-1]})
	stdinvy = crypto.Rebalance(cryptoParams, stdinvy)

	if debug && pid > 0 {
		SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{stdinvx}, nsnps, -1, ast.general.CachePath("stdinvx.txt"))                  // 1 / sqrt(sxx - (sx*sx/n))
		SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{crypto.CipherVector{stdinvy}}, 1, -1, ast.general.CachePath("stdinvy.txt")) // 1 / sqrt(syy - (sy*sy/n))
	}

	return stdinvx, stdinvy
}
