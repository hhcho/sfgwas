package gwas

import (
	"fmt"
	"github.com/hhcho/sfgwas/mpc"
	"math"
	"runtime"
	"strconv"
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

// Multiplies an encrypted matrix with a block of genotypes, optionally squares the result (square = true)
func (ast *AssocTest) GenoBlockMult(b int, mat crypto.CipherMatrix, square bool) (matOut crypto.CipherMatrix, dosageSum, dosageSqSum []float64, filtOut []bool) {
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

	// for logistic this method is used in multiple multiplications
	if !ast.general.config.UseLogistic && fileExists(multFile) && fileExists(dosFile) && fileExists(dos2File) && fileExists(filtFile) {

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

						var mult crypto.CipherMatrix
						var sum, sqSum []float64
						mult, _, _ = MatMult4Stream(cryptoParams, mat, X, 5, false, square, nprocsPerBlock)
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
			matOut, dosageSum, dosageSqSum = MatMult4Stream(cryptoParams, mat, XBlock, 5, true, square, 0)

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
	XBlocks := ast.general.genoBlocks

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
		// TODO include pid == 0 case in main workflow instead of separated

		if ast.general.config.UseLogistic {
			nbrCovTot := gwasParams.NumCov() + gwasParams.numPCs + 1
			// Party 0 is helping for logistic regression
			empty := mpc_core.InitRMat(ast.general.mpcObj[0].GetRType().Zero(), nbrCovTot, nbrCovTot)
			log.LLvl1("starts inverse")

			for i := 0; i < ast.general.config.Epochs+1; i++ { //+1 is for the last inversion done in the statistical test
				log.LLvl1("starts inverse new, iter ", i)
				_ = ast.general.mpcObj[0].MatrixInverseSqrtSVD(empty)
			}

			// helps computing operations in SS
			nbrOfBlocks := 0
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				newNbrOfBlocks := mpcObj.Network.ReceiveInt(party)
				if party != 1 {
					if newNbrOfBlocks != nbrOfBlocks {
						log.Fatal("Invalid number of SNPs among parties")
					}
				}
				nbrOfBlocks = newNbrOfBlocks

			}
			for block := 0; block < nbrOfBlocks; block++ {
				nbrOfSNPs := 0
				for party := 1; party <= ast.general.config.NumMainParties; party++ {
					newNbrOfSNPs := mpcObj.Network.ReceiveInt(party)
					if party != 1 {
						if newNbrOfSNPs != nbrOfSNPs {
							log.Fatal("Invalid number of SNPs among parties")
						}
					}
					nbrOfSNPs = newNbrOfSNPs

				}

				nbrOfSNPbatches := int(math.Ceil(float64(nbrOfSNPs) / float64(cryptoParams.GetSlots())))
				log.LLvl1("Number of SNPs: ", nbrOfSNPs, nbrOfBlocks, " Number of SNP batches: ", nbrOfSNPbatches)

				emptyVec := mpc_core.InitRVec(ast.general.mpcObj[0].GetRType().Zero(), nbrOfSNPs)
				_, _ = ast.general.mpcObj[0].SqrtAndSqrtInverse(emptyVec, true)

				_ = ast.general.mpcObj.SSMultElemVec(emptyVec, emptyVec)
				ast.general.mpcObj[0].TruncVec(emptyVec, ast.general.mpcObj[0].GetDataBits(), ast.general.mpcObj[0].GetFracBits())
				_ = ast.general.mpcObj.RevealSymVec(emptyVec)
			}
		} else {
			numCtx = mpcObj.Network.ReceiveInt(mpcObj.GetHubPid())
			nsnps = mpcObj.Network.ReceiveInt(mpcObj.GetHubPid())

			varx = crypto.CZeros(cryptoParams, numCtx)
			vary = crypto.CZeros(cryptoParams, 1)
		}

	} else { // pid > 0

		if ast.general.config.UseLogistic {
			approxInterval := crypto.IntervalApprox{A: ast.general.config.A, B: ast.general.config.B, Degree: ast.general.config.Degree}
			yClear := crypto.DecodeFloatVector(cryptoParams, y)

			timeTrain := time.Now()
			// compute covariates weights
			weights := crypto.CZeros(cryptoParams, ncov)
			if fileExists(ast.general.CachePath("finalWeights_party"+strconv.Itoa(pid)+".txt")) && debug {
				weights = LoadCacheFromFile(cryptoParams, ast.general.CachePath("finalWeights_party"+strconv.Itoa(pid)+".txt"))[0]
				// Scale by 1/sqrt(n) as the result of the QR is scaled UP by sqrt(n)
				log.LLvl1("Scaling Q by 1/sqrt(n)")
				Q = crypto.CMultConstMat(cryptoParams, Q, 1.0/math.Sqrt(float64(nrowsTotal)), false)
			} else {
				log.LLvl1("Scaling Q by 1/sqrt(n)")
				Q = crypto.CMultConstMat(cryptoParams, Q, 1.0/math.Sqrt(float64(nrowsTotal)), false)
				Q = crypto.CMatRescale(cryptoParams, Q)
				// Newton's method
				weights = ast.LrNewtonBasedCovOnly(Q, yClear, ast.general.config.Epochs, ast.general.config.LocalNumThreads,
					nrowsAll, approxInterval, 0.0)
			}
			log.LLvl1("Saving covariates-only weights" + strconv.Itoa(pid) + ".txt")
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				SaveMatrixToFile(ast.general.cps, ast.general.mpcObj[0], crypto.CipherMatrix{weights}, ncov,
					party, ast.general.CachePath("finalWeights_party"+strconv.Itoa(party)+".txt"))
			}
			log.LLvl1("Time for gradient descent and intercept: ", time.Since(timeTrain))
			ast.general.mpcObj.GetNetworks().PrintNetworkLog()

			// score test based on null model approach
			timeAll := time.Now()
			numInds := ast.general.config.NumInds[pid]
			tStatAllSNPs, outFilt := ast.gWY(Q, XBlocks, yClear, numInds, ast.general.config.LocalNumThreads, weights, approxInterval, nrowsAll)
			log.LLvl1("Time for score test part: ", time.Since(timeAll))

			log.LLvl1("LOGISTIC DONE")
			ast.general.mpcObj.GetNetworks().PrintNetworkLog()

			return tStatAllSNPs, outFilt

		} else {
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
					concatOut, dosageSum, dosageSqSum, filt := ast.GenoBlockMult(b, concat, false)
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
	}

	if !ast.general.config.UseLogistic {
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

func (ast *AssocTest) LrNewtonBasedCovOnly(C crypto.CipherMatrix, y []float64, iter int,
	numThreads int, numInds []int, approx crypto.IntervalApprox, initWeight float64) crypto.CipherVector {
	// init weights
	initClear := make([]float64, ast.general.cps.GetSlots())
	for i := range initClear {
		initClear[i] = initWeight
	}
	weights, _ := crypto.EncryptFloatVector(ast.general.cps, initClear)

	log.LLvl1("Scaling up the covariates matrix")
	totalInds := 0
	for _, v := range numInds {
		totalInds += v
	}
	CscaledUp := crypto.CMultConstMat(ast.general.cps, C, math.Sqrt(float64(totalInds)), false) // scale up to compute Wz
	CscaledUp = crypto.CMatRescale(ast.general.cps, CscaledUp)

	for i := 0; i < iter; i++ {
		grad, invHess, _, _, _, _ := ast.computeGradAndInvHessian(C, CscaledUp, y, weights, numThreads, approx, numInds, i == 0, true, i)
		step := CMultMatInnerProdVector(ast.general.cps, invHess, grad, len(C), numThreads)
		if ast.general.config.Debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				mpc.SaveMatrixToFileWithPrintIndex(ast.general.cps, ast.general.mpcObj[0], crypto.CipherMatrix{step}, len(C),
					party, ast.general.CachePath("step_"+strconv.Itoa(i)+".txt"), true, 0, len(C))
			}
		}
		weights = crypto.CAdd(ast.general.cps, weights, step)
		if ast.general.config.Debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				mpc.SaveMatrixToFileWithPrintIndex(ast.general.cps, ast.general.mpcObj[0], crypto.CipherMatrix{weights}, len(C),
					party, ast.general.CachePath("weights_"+strconv.Itoa(i)+".txt"), true, 0, len(C))
			}
		}
	}
	return weights
}

func (ast *AssocTest) computeGradAndInvHessian(C crypto.CipherMatrix, CScaledUp crypto.CipherMatrix, y []float64, weights crypto.CipherVector, numThreads int,
	approx crypto.IntervalApprox, nrowsAll []int, skipIntercept, computeGrad bool, iter int) (crypto.CipherVector,
	crypto.CipherMatrix, crypto.CipherMatrix, crypto.CipherMatrix, crypto.CipherVector, crypto.CipherVector) {

	totalnbrInds := 0
	for _, v := range nrowsAll {
		totalnbrInds += v
	}

	debug := ast.general.config.Debug
	mutex := sync.Mutex{}
	cryptoParams := ast.general.cps

	timepHatYtilde := time.Now()

	log.LLvl1("Computing u") // u in manuscript, phat in code
	maxInd := 0
	for m := 0; m < len(ast.general.config.NumInds); m++ {
		if ast.general.config.NumInds[m] > maxInd {
			maxInd = ast.general.config.NumInds[m]
		}
	}

	var covIntercept crypto.CipherVector
	if skipIntercept {
		zeros := make([]float64, len(y))
		covIntercept, _ = crypto.EncryptFloatVector(cryptoParams, zeros)
	} else {
		covIntercept = CMultMatColTimesColToCol(ast.general.cps, C, crypto.CipherMatrix{weights},
			ast.general.config.NumInds[ast.general.mpcObj[0].GetPid()], numThreads)[0]
	}
	if debug {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{covIntercept}, nrowsAll[party],
				party, ast.general.CachePath("covIntercept"+strconv.Itoa(iter)+".txt"), true)
		}
	}

	covInterceptPadded := make(crypto.CipherVector, int(math.Ceil(float64(maxInd)/float64(ast.general.cps.GetSlots()))))
	for f := range covInterceptPadded {
		if f < len(covIntercept) {
			covInterceptPadded[f] = covIntercept[f]
		} else {
			covInterceptPadded[f] = covIntercept[0].CopyNew().Ciphertext()
		}

	}
	if debug {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{covInterceptPadded}, nrowsAll[party],
				party, ast.general.CachePath("pHat"+strconv.Itoa(iter)+".txt"), true)
		}
	}
	pHatPadded := mpc.CSigmoidApprox(ast.general.mpcObj[0].GetPid(), ast.general.mpcObj[0].Network, cryptoParams, covInterceptPadded,
		approx, &mutex) // n x 1

	pHat := pHatPadded[:len(covIntercept)]
	pHat = ast.general.mpcObj[0].Network.BootstrapMatAll(cryptoParams, crypto.CipherMatrix{pHat})[0]

	if debug {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{pHat}, nrowsAll[party],
				party, ast.general.CachePath("pHat"+strconv.Itoa(iter)+".txt"), true)
		}
	}

	log.LLvl1("Computing yTilde")
	yEncoded, _ := crypto.EncodeFloatVector(cryptoParams, y)  // n x 1
	yTilde := crypto.CPSubOther(cryptoParams, yEncoded, pHat) // n x 1

	if debug {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			mpc.SaveMatrixToFileWithPrint(ast.general.cps, ast.general.mpcObj[0], crypto.CipherMatrix{yTilde}, nrowsAll[party],
				party, ast.general.CachePath("yTilde"+strconv.Itoa(iter)+".txt"), true)
		}
	}
	var grad crypto.CipherVector // r in manuscript, grad in code
	if computeGrad {
		grad = CMultMatInnerProdVector(cryptoParams, C, yTilde, ast.general.config.NumInds[ast.general.mpcObj[0].GetPid()], numThreads)
		grad = ast.general.mpcObj[0].Network.AggregateCVec(cryptoParams, grad)
		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				SaveMatrixToFile(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{grad}, len(C),
					party, ast.general.CachePath("grad"+strconv.Itoa(iter)+".txt"))
			}
		}
	}

	log.LLvl1("Computing w") // Diagonal of R in manuscript, w in code
	pHatSquare := crypto.CMult(cryptoParams, pHat, pHat)

	w := crypto.CSub(cryptoParams, pHat, pHatSquare)
	w = ast.general.mpcObj[0].Network.BootstrapMatAll(cryptoParams, crypto.CipherMatrix{w})[0]
	if debug {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			SaveMatrixToFile(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{w}, nrowsAll[party],
				party, ast.general.CachePath("w"+strconv.Itoa(iter)+".txt"))
		}
	}

	log.LLvl1("Time to compute u and yTilde: ", time.Since(timepHatYtilde))

	timeWz := time.Now()

	log.LLvl1("Computing V (using CScaledUp)") // V in manuscript, Wz in code
	Wz := make(crypto.CipherMatrix, len(C))    // column-encrypted, n x c

	for i := 0; i < len(C); i++ {
		Wz[i] = crypto.CMult(cryptoParams, w, CScaledUp[i]) // Wz is scaled up
	}
	if debug {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			SaveMatrixToFile(cryptoParams, ast.general.mpcObj[0], Wz, nrowsAll[party],
				party, ast.general.CachePath("Wz"+strconv.Itoa(iter)+".txt"))
		}
	}
	log.LLvl1("Time to compute Wz: ", time.Since(timeWz))
	timeZTwZInv := time.Now()

	var ZTwZInv crypto.CipherMatrix
	log.LLvl1("Computing W ")                                  // W in manuscript, ZTwZInv in code
	ZTwZ := CMultMatInnerProd(cryptoParams, C, Wz, numThreads) // C = ZT row-wise encrypted, result is row-wise, #cov x #cov
	// compute (ZTwZ)^(-1)
	ZTwZ = ast.general.mpcObj[0].Network.AggregateCMat(cryptoParams, ZTwZ)

	ZTwZ = ast.general.mpcObj[0].Network.CollectiveBootstrapMat(cryptoParams, ZTwZ, -1)

	// TODO include scaling factor as a parameter
	ZTwZ = crypto.CMultConstMat(cryptoParams, ZTwZ, 1.0/(float64(totalnbrInds)/ast.general.config.InverseMatScale), false) // unscale from scale Wz
	ZTwZ = crypto.CMatRescale(cryptoParams, ZTwZ)

	if debug {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], ZTwZ, len(C),
				party, ast.general.CachePath("ZTwZ"+strconv.Itoa(iter)+".txt"), true)
		}
	}

	log.LLvl1("Computing inverse of W using SS")

	ZTwZSS := ast.general.mpcObj[0].CMatToSS(cryptoParams, ast.general.mpcObj[0].GetRType(), ZTwZ, -1, len(ZTwZ), 1, len(C))

	log.LLvl1("converted W to SS, getting BT such as BTB = inverse(ZTwZ)")
	BSS := ast.general.mpcObj[0].MatrixInverseSqrtSVD(ZTwZSS)
	BT := ast.general.mpcObj[0].SSToCMat(cryptoParams, BSS.Transpose())

	// scale back such that BTB ends up unscaled
	BT = crypto.CMultConstMat(cryptoParams, BT, math.Sqrt((ast.general.config.InverseMatScale*2)/math.Sqrt(float64(totalnbrInds))), false)
	BT = crypto.CMatRescale(cryptoParams, BT)

	// TODO check if can be removed
	maskClear := make([]float64, ast.general.cps.GetSlots())
	for i := 0; i < len(C); i++ {
		maskClear[i] = 1
	}
	maskEnc, _ := crypto.EncodeFloatVectorWithScale(cryptoParams, maskClear, float64(ast.general.cps.Params.Qi()[BT[0][0].Level()]))
	for i := range BT {
		BT[i] = crypto.CPMult(cryptoParams, BT[i], maskEnc)
	}

	ZTwZInv = CMultMatInnerProd(cryptoParams, BT, crypto.CopyEncryptedMatrix(BT), numThreads)

	if debug {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], ZTwZInv, len(C),
				party, ast.general.CachePath("ZTwZInv"+strconv.Itoa(iter)+".txt"), true)
		}
	}
	log.LLvl1("Time to compute ZTwZInv: ", time.Since(timeZTwZInv))

	if debug {
		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], BT, len(C),
				party, ast.general.CachePath("ZTwZInvNew"+strconv.Itoa(iter)+".txt"), true)
		}
	}
	log.LLvl1("Time to compute ZTwZInvNew: ", time.Since(timeZTwZInv))

	return grad, ZTwZInv, BT, Wz, yTilde, w
}

// gWY computes the score test for logistic regression
func (ast *AssocTest) gWY(C crypto.CipherMatrix, XBlocks []*GenoFileStream, y []float64, numInds int, numThreads int,
	weights crypto.CipherVector, approx crypto.IntervalApprox, nrowsAll []int) (crypto.CipherVector, []bool) {

	totalInds := 0
	for _, v := range nrowsAll {
		totalInds += v
	}

	debug := ast.general.config.Debug
	cryptoParams := ast.general.cps

	CscaledUp := crypto.CMultConstMat(ast.general.cps, C, math.Sqrt(float64(totalInds)), false) // scale up to compute Wz
	CscaledUp = crypto.CMatRescale(ast.general.cps, CscaledUp)

	// common workflow as in newton's method
	_, ZTwZInv, BT, Wz, yTilde, w := ast.computeGradAndInvHessian(C, CscaledUp, y, weights, numThreads,
		approx, nrowsAll, false, false, -1)

	// computes the different elements of the statistical test
	timeWzZTwZInv := time.Now()

	log.LLvl1("Computing U (WzZTwZInv)")
	var WzZTwZInv crypto.CipherMatrix
	if fileExists(ast.general.CachePath("WzZTwZInv.txt")) && debug {
		WzZTwZInv = LoadCacheFromFile(cryptoParams, ast.general.CachePath("WzZTwZInv.txt"))
	} else {
		WzZTwZInv = CMultMatColTimesRowToCol(cryptoParams, Wz, ZTwZInv, numInds, len(C), numThreads)
		WzZTwZInv = ast.general.mpcObj[0].Network.BootstrapMatAll(cryptoParams, WzZTwZInv)

		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				SaveMatrixToFile(cryptoParams, ast.general.mpcObj[0], WzZTwZInv, nrowsAll[party],
					party, ast.general.CachePath("WzZTwZInv.txt"))
			}
		}
	}
	log.LLvl1("Time to compute WzZTwZInv: ", time.Since(timeWzZTwZInv))

	timeWzBT := time.Now()
	log.LLvl1("Computing E (WzBT)")
	var WzBT crypto.CipherMatrix
	if fileExists(ast.general.CachePath("WzBT.txt")) && debug {
		WzBT = LoadCacheFromFile(cryptoParams, ast.general.CachePath("WzBT.txt"))
	} else {
		WzBT = CMultMatColTimesRowToCol(cryptoParams, Wz, BT, numInds, len(C), numThreads)
		WzBT = ast.general.mpcObj[0].Network.BootstrapMatAll(cryptoParams, WzBT)

		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				SaveMatrixToFile(cryptoParams, ast.general.mpcObj[0], WzBT, nrowsAll[party],
					party, ast.general.CachePath("WzBT.txt"))
			}
		}
	}
	log.LLvl1("Time to compute E (WzBT): ", time.Since(timeWzBT))
	timeWzZTwZInvZTy := time.Now()

	log.LLvl1("Computing o (WzZTwZInvZTy)")
	var WzZTwZInvZTy crypto.CipherVector
	if fileExists(ast.general.CachePath("WzZTwZInvZTy.txt")) && debug {
		WzZTwZInvZTy = LoadCacheFromFile(cryptoParams, ast.general.CachePath("WzZTwZInvZTy.txt"))[0]
	} else {
		log.LLvl1("Computing d (ZTy), a part of o")
		mask := make([]float64, ast.general.cps.GetSlots()*len(C[0]))
		for i := 0; i < numInds; i++ {
			mask[i] = 1
		}
		maskEncoded, _ := crypto.EncodeFloatVector(ast.general.cps, mask)
		yTilde = crypto.CPMult(cryptoParams, yTilde, maskEncoded)

		yTilde = ast.general.mpcObj[0].Network.BootstrapMatAll(cryptoParams, crypto.CipherMatrix{yTilde})[0]
		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{yTilde}, numInds,
					party, ast.general.CachePath("yTilde_befZTy.txt"), true)
				mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], CscaledUp, numInds,
					party, ast.general.CachePath("CScaledup_befZTy.txt"), true)
			}
		}

		ZTy := CMultMatInnerProdVector(cryptoParams, CscaledUp, yTilde, numInds, numThreads) // ZTy is a vector of #cov elements
		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{ZTy}, len(C),
					party, ast.general.CachePath("ZTy_local_p_"+strconv.Itoa(party)+".txt"), true)
			}
		}
		ZTy = ast.general.mpcObj[0].Network.AggregateCVec(cryptoParams, ZTy)
		ZTy = ast.general.mpcObj[0].Network.CollectiveBootstrapVec(cryptoParams, ZTy, -1)

		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{ZTy}, len(C),
					party, ast.general.CachePath("ZTy.txt"), false)
			}
		}
		// need to do WzZTwZInv * ZTy
		WzZTwZInvZTy = CMultMatColTimesColToCol(cryptoParams, WzZTwZInv, crypto.CipherMatrix{ZTy}, numInds, numThreads)[0] // still scaled up, result is a vector of n elements

		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{WzZTwZInvZTy}, nrowsAll[party],
					party, ast.general.CachePath("WzZTwZInvZTyScaled.txt"), false)
			}
		}

		WzZTwZInvZTyScaleBackDown := crypto.CMultConstMat(cryptoParams, crypto.CipherMatrix{WzZTwZInvZTy}, (1.0 / float64(totalInds)), false)
		WzZTwZInvZTy = crypto.CMatRescale(cryptoParams, WzZTwZInvZTyScaleBackDown)[0]

		WzZTwZInvZTy = ast.general.mpcObj[0].Network.BootstrapMatAll(cryptoParams, crypto.CipherMatrix{WzZTwZInvZTy})[0]

		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{WzZTwZInvZTy}, nrowsAll[party],
					party, ast.general.CachePath("WzZTwZInvZTy.txt"), false)
			}
		}
	}

	log.LLvl1("Time to compute WzZTwZInvZTy: ", time.Since(timeWzZTwZInvZTy))

	ast.general.mpcObj.GetNetworks().PrintNetworkLog()

	// for all SNPs (send to p0 for its help)
	ast.general.mpcObj[0].Network.SendInt(len(XBlocks), 0)
	resultAllSnps := make(crypto.CipherVector, 0)
	outFilt := make([]bool, 0)

	for block := 0; block < len(XBlocks); block++ {
		timeBlock := time.Now()
		log.LLvl1("BLOCK: ", block)
		X := XBlocks[block]

		nbrOfSNPs := 0
		if !ast.general.IsPgen() {
			nbrOfSNPs = len(X.readRow())
		} else {
			blockSize := ast.general.genoBlockSizes[block]

			shift := uint64(0)
			for i := 0; i < block; i++ {
				shift += uint64(ast.general.genoBlockSizes[i])
			}

			if ast.general.gwasParams.snpFilt == nil { // no QC
				nbrOfSNPs = blockSize
			} else {
				nbrOfSNPs = SumBool(ast.general.gwasParams.snpFilt[shift : shift+uint64(blockSize)])
			}
		}

		ast.general.mpcObj[0].Network.SendInt(nbrOfSNPs, 0)
		nbrOfSNPciphers := int(math.Ceil(float64(nbrOfSNPs) / float64(cryptoParams.GetSlots())))
		log.LLvl1("Number of SNPs: ", nbrOfSNPs, " Number of SNP batches: ", nbrOfSNPciphers)

		timegTWzZTwZInvzTWg := time.Now()
		log.LLvl1("Computing b (gTWzZTwZInvzTWg)")
		zTWBTg := crypto.CZeroMat(cryptoParams, 1, len(C))
		gTWzZTwZInvzTWg := crypto.CZeros(cryptoParams, nbrOfSNPciphers)

		if fileExists(ast.general.CachePath("gTWzZTwZInvzTWg"+strconv.Itoa(block)+".txt")) && debug {
			gTWzZTwZInvzTWg = LoadCacheFromFile(cryptoParams, ast.general.CachePath("gTWzZTwZInvzTWg"+strconv.Itoa(block)+".txt"))[0]
		} else {
			// compute H
			zTWBTg, _, _, _ = ast.GenoBlockMult(block, WzBT, false) // WzBT is column encrypted so here we take its transpose
			if debug {
				for party := 1; party <= ast.general.config.NumMainParties; party++ {
					SaveMatrixToFile(cryptoParams, ast.general.mpcObj[0], zTWBTg, nbrOfSNPs,
						party, ast.general.CachePath("zTWBTg"+strconv.Itoa(block)+".txt"))
				}
			}
			zTWBTgAggr := ast.general.mpcObj[0].Network.AggregateCMat(cryptoParams, zTWBTg)
			zTWBTgAggr = ast.general.mpcObj[0].Network.CollectiveBootstrapMat(cryptoParams, zTWBTgAggr, -1)

			if debug {
				for party := 1; party <= ast.general.config.NumMainParties; party++ {
					SaveMatrixToFile(cryptoParams, ast.general.mpcObj[0], zTWBTgAggr, nbrOfSNPs,
						party, ast.general.CachePath("zTWBTgAggr"+strconv.Itoa(block)+".txt"))
				}
			}
			for i := 0; i < len(zTWBTg); i++ {
				toAdd := crypto.CMult(cryptoParams, crypto.CopyEncryptedVector(zTWBTgAggr[i]), crypto.CopyEncryptedVector(zTWBTgAggr[i]))
				gTWzZTwZInvzTWg = crypto.CAdd(cryptoParams, gTWzZTwZInvzTWg, toAdd)
			}
			gTWzZTwZInvzTWg = crypto.CMultConst(cryptoParams, gTWzZTwZInvzTWg, 1.0/float64(totalInds), false) // unscale from Wz
			gTWzZTwZInvzTWg = crypto.CRescale(cryptoParams, gTWzZTwZInvzTWg)
		}
		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				SaveMatrixToFile(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{gTWzZTwZInvzTWg}, nbrOfSNPs,
					party, ast.general.CachePath("gTWzZTwZInvzTWg"+strconv.Itoa(block)+".txt"))
			}
		}
		log.LLvl1("Time to compute b (gTWzZTwZInvzTWg): ", time.Since(timegTWzZTwZInvzTWg))

		timegTWg := time.Now()
		var gTWg crypto.CipherMatrix
		log.LLvl1("Computing x (gTWg)")
		if fileExists(ast.general.CachePath("gTWg"+strconv.Itoa(block)+".txt")) && debug {
			gTWg = LoadCacheFromFile(cryptoParams, ast.general.CachePath("gTWg"+strconv.Itoa(block)+".txt"))
		} else {
			gTWg, _, _, _ = ast.GenoBlockMult(block, crypto.CipherMatrix{w}, true)

			if debug {
				for party := 1; party <= ast.general.config.NumMainParties; party++ {
					mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], gTWg, nbrOfSNPs,
						party, ast.general.CachePath("gTWg"+strconv.Itoa(block)+".txt"), true)
				}
			}
		}
		log.LLvl1("Time to compute x (gTWg): ", time.Since(timegTWg))

		gTWgAggr := ast.general.mpcObj[0].Network.AggregateCMat(cryptoParams, gTWg)
		gTildesTgTilde := crypto.CSub(cryptoParams, gTWgAggr[0], gTWzZTwZInvzTWg)
		// ensures that all parties have the input ciphertext for the following step
		gTildesTgTilde = ast.general.mpcObj[0].Network.CollectiveBootstrapMat(cryptoParams, crypto.CipherMatrix{gTildesTgTilde}, 1)[0]

		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{gTildesTgTilde}, nbrOfSNPs,
					party, ast.general.CachePath("gTildesTgTildeAggr_"+strconv.Itoa(block)+".txt"), true)
			}
		}

		timegTyTilde := time.Now()
		var gTyTilde crypto.CipherMatrix
		log.LLvl1("Computing j (gTyTilde)")
		if fileExists(ast.general.CachePath("gTyTilde"+strconv.Itoa(ast.general.mpcObj[0].GetPid())+"_"+strconv.Itoa(block)+".txt")) && debug {
			gTyTilde = LoadCacheFromFile(cryptoParams, ast.general.CachePath("gTyTilde"+strconv.Itoa(ast.general.mpcObj[0].GetPid())+"_"+strconv.Itoa(block)+".txt"))
		} else {
			gTyTilde, _, _, _ = ast.GenoBlockMult(block, crypto.CipherMatrix{yTilde}, false)
			log.LLvl1("Time to compute j (gTyTilde): ", time.Since(timegTyTilde))
			if debug {
				for party := 1; party <= ast.general.config.NumMainParties; party++ {
					mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], gTyTilde, nbrOfSNPs,
						party, ast.general.CachePath("gTyTilde"+strconv.Itoa(party)+"_"+strconv.Itoa(block)+".txt"), false)
				}
			}
		}

		var filtOutBlock []bool
		timegTWzZTwZInvZTy := time.Now()
		var gTWzZTwZInvZTy crypto.CipherMatrix
		log.LLvl1("Computing s (gTWzZTwZInvZTy)")
		if fileExists(ast.general.CachePath("gTWzZTwZInvZTy"+strconv.Itoa(ast.general.mpcObj[0].GetPid())+"_"+strconv.Itoa(block)+".txt")) && debug {
			gTWzZTwZInvZTy = LoadCacheFromFile(cryptoParams, ast.general.CachePath("gTWzZTwZInvZTy"+strconv.Itoa(ast.general.mpcObj[0].GetPid())+"_"+strconv.Itoa(block)+".txt"))
		} else {
			// filtout comes from linear workflow for comptability
			gTWzZTwZInvZTy, _, _, filtOutBlock = ast.GenoBlockMult(block, crypto.CipherMatrix{WzZTwZInvZTy}, false)
			log.LLvl1("Time to compute s (gTWzZTwZInvZTy): ", time.Since(timegTWzZTwZInvZTy))
			if debug {
				for party := 1; party <= ast.general.config.NumMainParties; party++ {
					mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], gTWzZTwZInvZTy, nbrOfSNPs,
						party, ast.general.CachePath("gTWzZTwZInvZTy"+strconv.Itoa(party)+"_"+strconv.Itoa(block)+".txt"), false)
				}
			}
		}

		log.LLvl1("Bootstrapping and Aggregating --> f (gTildesTyTildeAggr)")
		gTyTilde = ast.general.mpcObj[0].Network.BootstrapMatAll(cryptoParams, gTyTilde)
		gTWzZTwZInvZTy = ast.general.mpcObj[0].Network.BootstrapMatAll(cryptoParams, gTWzZTwZInvZTy)

		gTildesTyTilde := crypto.CSub(cryptoParams, gTyTilde[0], gTWzZTwZInvZTy[0])
		gTildesTyTilde = crypto.CRescale(cryptoParams, gTildesTyTilde)
		gTildesTyTildeAggr := ast.general.mpcObj[0].Network.AggregateCVec(cryptoParams, gTildesTyTilde)

		if debug {
			for party := 1; party <= ast.general.config.NumMainParties; party++ {
				mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{gTildesTyTildeAggr}, nbrOfSNPs,
					party, ast.general.CachePath("gTildesTyTildeAggr"+strconv.Itoa(block)+".txt"), false)
			}
		}

		timeFinalStep := time.Now()

		gTildesTyTildeAggrSS := mpc_core.InitRVec(ast.general.mpcObj[0].GetRType().Zero(), cryptoParams.GetSlots()*nbrOfSNPciphers)
		gTildesTgTildeAggrSS := mpc_core.InitRVec(ast.general.mpcObj[0].GetRType().Zero(), cryptoParams.GetSlots()*nbrOfSNPciphers)

		log.LLvl1("Transforming numerator f to SS")
		gTildesTyTildeAggrSS = ast.general.mpcObj[0].CVecToSS(cryptoParams, ast.general.mpcObj[0].GetRType(), gTildesTyTildeAggr, -1, nbrOfSNPciphers, nbrOfSNPs)
		log.LLvl1("Transforming denominator j to SS")
		gTildesTgTildeAggrSS = ast.general.mpcObj[0].CVecToSS(cryptoParams, ast.general.mpcObj[0].GetRType(), gTildesTgTilde, -1, nbrOfSNPciphers, nbrOfSNPs)

		log.LLvl1("Sqrt Inverse of denominator j in SS")
		_, gTildesTgTildeAggrSSInvSqrt := ast.general.mpcObj[0].SqrtAndSqrtInverse(gTildesTgTildeAggrSS, true)

		divOut := ast.general.mpcObj.SSMultElemVec(gTildesTyTildeAggrSS, gTildesTgTildeAggrSSInvSqrt)
		divOut = ast.general.mpcObj[0].TruncVec(divOut, ast.general.mpcObj[0].GetDataBits(), ast.general.mpcObj[0].GetFracBits())

		outCipher := ast.general.mpcObj[0].SSToCVec(cryptoParams, divOut)
		log.LLvl1("time for final step: ", time.Since(timeFinalStep))
		log.LLvl1("time for block: ", time.Since(timeBlock))

		for party := 1; party <= ast.general.config.NumMainParties; party++ {
			mpc.SaveMatrixToFileWithPrint(cryptoParams, ast.general.mpcObj[0], crypto.CipherMatrix{outCipher}, nbrOfSNPs,
				party, ast.general.CachePath("outCipher_"+strconv.Itoa(block)+".txt"), true)
		}
		resultAllSnps = append(resultAllSnps, outCipher...)
		outFilt = append(outFilt, filtOutBlock...)
		ast.general.mpcObj.GetNetworks().PrintNetworkLog()
	}
	return resultAllSnps, outFilt
}
