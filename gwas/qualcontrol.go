package gwas

import (
	// "math"

	"fmt"
	"math"
	"runtime"
	"time"

	"go.dedis.ch/onet/v3/log"

	mpc_core "github.com/hhcho/mpc-core"
	"gonum.org/v1/gonum/mat"
)

type QC struct {
	general *ProtocolInfo

	filtNumInds []int //total number of individuals that pass quality control filters
	filtNumSnps int   //total number of SNPs that pass quality control filters

	filterParams *FilterParams
}

func (g *ProtocolInfo) InitQC(filterParams *FilterParams) QC {
	return QC{
		general:      g,
		filtNumInds:  make([]int, g.mpcObj[0].GetNParty()),
		filtNumSnps:  0,
		filterParams: filterParams,
	}
}

//IndividualMissAndHetFilters filters individuals based on missing rate and heterozygosity filter
func (qc *QC) IndividualMissAndHetFilters() []bool {
	if qc.general.mpcObj[0].GetPid() == 0 {
		return make([]bool, 1)
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Computing local individual filters (missing rate and heterozygosity)")

	numInds := qc.general.genoBlocks[0].NumRows()
	miss := make([]int, numInds)
	het := make([]int, numInds)
	for _, genoFs := range qc.general.genoBlocks { // TODO: parallelize
		for row, idx := genoFs.NextRow(), 0; row != nil; row, idx = genoFs.NextRow(), idx+1 {
			for i := range row {
				if row[i] < 0 {
					miss[idx]++
				}
				if row[i] == 1 {
					het[idx]++
				}
			}
		}

		genoFs.Reset()
	}

	numSnps := qc.filtNumSnps
	ikeep := make([]bool, numInds)
	for i := range ikeep {
		missRate := float64(miss[i]) / float64(numSnps)
		hetRate := float64(het[i]) / float64(int(numSnps)-miss[i])

		//if i < 100 {
		//	log.LLvl1(time.Now().Format(time.RFC3339), "missRate(%f) missBound(%f) hetRate(%f) hetBound(%f,%f)\n", missRate, qc.filterParams.IndMissBound, hetRate, qc.filterParams.HetLowerBound, qc.filterParams.HetUpperBound)
		//}

		if missRate < qc.filterParams.IndMissBound &&
			hetRate < qc.filterParams.HetUpperBound &&
			hetRate > qc.filterParams.HetLowerBound {
			ikeep[i] = true
		} else {
			ikeep[i] = false
		}
	}

	return ikeep
}

// ac: allele counts (0, 1)
// gc: genotype counts (0, 1, 2)
// miss: missing value counts
func (qc *QC) SNPFilterWithPrecomputedStats(ac, gc [][]uint32, miss []uint32, useCache bool) []bool {
	mpcPar := qc.general.mpcObj
	pid := mpcPar[0].GetPid()
	numSnpWindow := len(miss)
	maxWindowSize := 10000000

	if numSnpWindow > maxWindowSize { // Work on maxWindowSize snps at a time
		out := make([]bool, numSnpWindow)
		for start, end, batchIndex := 0, maxWindowSize, 0; start < numSnpWindow; start, end, batchIndex = start+maxWindowSize, end+maxWindowSize, batchIndex+1 {
			if end > numSnpWindow {
				end = numSnpWindow
			}

			// Take a subset of input data
			acSub := make([][]uint32, len(ac))
			for i := range acSub {
				acSub[i] = ac[i][start:end]
			}
			gcSub := make([][]uint32, len(gc))
			for i := range gcSub {
				gcSub[i] = gc[i][start:end]
			}
			missSub := miss[start:end]

			// Run QC on the subset
			startTime := time.Now()
			log.LLvl1(startTime, fmt.Sprintf("Variant QC: started processing %d-%d / %d", start, end, numSnpWindow))
			outSub := qc.SNPFilterWithPrecomputedStats(acSub, gcSub, missSub, useCache)
			runtime.GC() // Clean up memory
			log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Variant QC: finished processing %d-%d / %d,", start, end, numSnpWindow), time.Since(startTime))

			if pid > 0 {
				cacheFile := qc.general.CachePath(fmt.Sprintf("gkeep.%d.bin", batchIndex))
				writeFilterToFile(cacheFile, outSub, true)
				log.LLvl1(time.Now().Format(time.RFC3339), "QC filter wrote to cache:", cacheFile)
			}

			copy(out[start:end], outSub)
		}
		return out
	}

	rtype := mpcPar[0].GetRType().Zero()
	useBoolean := mpcPar[0].GetBooleanShareFlag()

	gwasParams := qc.general.gwasParams
	// numSnp := gwasParams.numSnps
	numSnp := numSnpWindow
	numInd := uint32(gwasParams.numInds[pid])

	log.LLvl1(time.Now().Format(time.RFC3339), "Computing SNP missingness filter")

	totalInds := Sum(gwasParams.NumInds())

	var jkeep []bool
	{ // Scope for missingness filter
		xCount := mpc_core.InitRVec(rtype, numSnp)
		if pid > 0 {
			for i := range xCount {
				xCount[i] = rtype.FromUint64(uint64(numInd - miss[i]))
			}
		}

		/* SNP missing rate */
		lb := int((1 - qc.filterParams.GenoMissBound) * float64(totalInds))

		start := time.Now()
		log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Secure comparison for SNP missingness (m = %d) ... ", len(xCount)))
		gmissFilt := mpcPar.NotLessThanPublic(xCount, rtype.FromInt(lb), useBoolean) // Secure comparison
		log.LLvl1(time.Now().Format(time.RFC3339), "done. ", time.Since(start))

		snpFilt := mpcPar.RevealSymVec(gmissFilt)

		jkeep = make([]bool, len(snpFilt))
		if pid > 0 {
			for i := range jkeep {
				jkeep[i] = snpFilt[i].Uint64() != 0
			}
		}
	}
	runtime.GC()

	var numSnpKeep int
	if pid > 0 {
		numSnpKeep = SumBool(jkeep)
		if pid == mpcPar[0].GetHubPid() {
			mpcPar[0].Network.SendInt(numSnpKeep, 0)
		}
	} else if pid == 0 {
		numSnpKeep = mpcPar[0].Network.ReceiveInt(mpcPar[0].GetHubPid())
	}
	log.LLvl1(time.Now().Format(time.RFC3339), "Number of SNPs remaining after missingness filter:", numSnpKeep)

	log.LLvl1(time.Now().Format(time.RFC3339), "Computing SNP minor allele frequency filter")

	fracBits := mpcPar[0].GetFracBits()
	mafBound := qc.filterParams.MafLowerBound
	hweBound := rtype.FromFloat64(qc.filterParams.HweUpperBound, fracBits)

	/* Minor allele frequency (MAF) */
	{ // Scope for MAF filter

		// Target: bound <= s/c <= 1 - bound
		// Compute as: (2*s-c)^2 <= c^2 * (2*bound - 1)^2
		xCount := mpc_core.InitRVec(rtype, numSnpKeep)
		xSum := mpc_core.InitRVec(rtype, numSnpKeep) // Total count of non-reference allele
		xSumMinusCount := mpc_core.InitRVec(rtype, numSnpKeep)

		if pid > 0 {
			index := 0
			for i := 0; i < numSnp; i++ {
				if jkeep[i] {
					xCount[index] = rtype.FromUint64(2 * uint64(numInd-miss[i]))
					xSum[index] = rtype.FromUint64(uint64(ac[1][i]))
					index++
				}
			}
			xSumMinusCount = xSum.Copy()
			xSumMinusCount.MulScalar(rtype.FromInt(2))
			xSumMinusCount.Sub(xCount)
		}

		xSumSq := mpcPar.SSSquareElemVec(xSumMinusCount)
		xCountSq := mpcPar.SSSquareElemVec(xCount)

		if pid > 0 {
			prec := 20
			bound := rtype.FromFloat64(math.Pow(2*mafBound-1.0, 2), prec)

			xCountSq.MulScalar(bound)
			xSumSq.MulScalar(rtype.FromUint64(uint64(1) << prec))

			xCountSq.Sub(xSumSq)
		}

		start := time.Now()
		log.LLvl1(time.Now().Format(time.RFC3339), "IsPositive for MAF filter ... ")
		mafFilt := mpcPar.IsPositive(xCountSq, useBoolean)
		log.LLvl1(time.Now().Format(time.RFC3339), "done.", time.Since(start))

		snpFilt := mpcPar.RevealSymVec(mafFilt)

		if pid > 0 {
			index := 0
			for i := range jkeep {
				if jkeep[i] {
					jkeep[i] = snpFilt[index].Uint64() != 0
					index += 1
				}
			}
		}
	}
	runtime.GC()

	if pid > 0 {
		numSnpKeep = SumBool(jkeep)
		if pid == mpcPar[0].GetHubPid() {
			mpcPar[0].Network.SendInt(numSnpKeep, 0)
		}
	} else if pid == 0 {
		numSnpKeep = mpcPar[0].Network.ReceiveInt(mpcPar[0].GetHubPid())
	}
	log.LLvl1(time.Now().Format(time.RFC3339), "Number of SNPs remaining after MAF filter:", numSnpKeep)

	log.LLvl1(time.Now().Format(time.RFC3339), "Computing SNP Hardy-Weinberg equilibrium filter")

	/* Hardy-Weinberg equlibrium (over control cohort only) */
	// TODO: Using all subjects for now; for continuous phenotypes
	{ // Scope for HWE filter
		xSumCtrl := mpc_core.InitRVec(rtype, numSnpKeep)            // alpha: dosage sum
		xCountCtrl := mpc_core.InitRVec(rtype, numSnpKeep)          // beta: 2 * num observed samples
		genoObservedCtrl := mpc_core.InitRMat(rtype, 3, numSnpKeep) // 2 * genotype count
		if pid > 0 {
			index := 0
			for j := 0; j < numSnp; j++ {
				if jkeep[j] {
					xSumCtrl[index] = rtype.FromUint64(uint64(ac[1][j]))
					xCountCtrl[index] = rtype.FromUint64(2 * uint64(numInd-miss[j]))
					for i := range genoObservedCtrl {
						genoObservedCtrl[i][index] = rtype.FromUint64(uint64(gc[i][j]))
					}
					index++
				}
			}
		}

		// Calculate expected genotype frequencies eAA, eAa, eaa
		// alpha = dosage sum over observed controls
		// beta = 2 * observed controls size
		// eaa = (beta - alpha)^2 / (2 * beta)
		// eAa = (2 * alpha * (beta - alpha)) / (2 * beta)
		// eAA = alpha^2 / (2 * beta)
		// Note: Denominator multiplied to the bound
		xCountCtrlConst := xCountCtrl.Copy()
		xCountCtrlConst.MulScalar(rtype.FromInt(2)) // 2 * beta
		xCountCtrl.Sub(xSumCtrl)                    // beta - alpha

		expected := make([]mpc_core.RVec, 3)           // AA, Aa, aa
		expected[2] = mpcPar.SSSquareElemVec(xSumCtrl) // alpha^2
		xSumCtrl.MulScalar(rtype.FromInt(2))
		expected[1] = mpcPar.SSMultElemVec(xSumCtrl, xCountCtrl) // 2 * alpha * (beta - alpha)
		expected[0] = mpcPar.SSSquareElemVec(xCountCtrl)         // (beta - alpha)^2

		// Compute chi-squared statistics
		// Sum { (genoObs * (2 * beta) - expected)^2 / expected }
		chiSq := mpc_core.InitRVec(rtype.Zero(), numSnpKeep)
		for i := range expected {
			tmp := mpcPar.SSMultElemVec(xCountCtrlConst, genoObservedCtrl[i])
			tmp.Sub(expected[i])
			tmp = mpcPar.SSSquareElemVec(tmp)

			start := time.Now()

			divOut := mpcPar.Divide(tmp, expected[i], useBoolean)
			chiSq.Add(divOut)

			log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Division round %d/%d complete, %s", i+1, len(expected), time.Since(start)))
		}

		xCountCtrlConst.MulScalar(hweBound)

		start := time.Now()
		log.LLvl1(time.Now().Format(time.RFC3339), "Secure comparison for HWE ... ")
		hweFilt := mpcPar.LessThan(chiSq, xCountCtrlConst, useBoolean)
		log.LLvl1(time.Now().Format(time.RFC3339), "done.", time.Since(start))

		snpFilt := mpcPar.RevealSymVec(hweFilt)

		if pid > 0 {
			index := 0
			for i := range jkeep {
				if jkeep[i] {
					jkeep[i] = snpFilt[index].Uint64() != 0
					index += 1
				}
			}
		}
	}
	runtime.GC()

	if pid > 0 {
		numSnpKeep = SumBool(jkeep)
		if pid == mpcPar[0].GetHubPid() {
			mpcPar[0].Network.SendInt(numSnpKeep, 0)
		}
	} else if pid == 0 {
		numSnpKeep = mpcPar[0].Network.ReceiveInt(mpcPar[0].GetHubPid())
	}
	log.LLvl1(time.Now().Format(time.RFC3339), "Number of SNPs remaining after HWE filter:", numSnpKeep)

	return jkeep
}

func (qc *QC) SNPMissFilter() []bool {

	log.LLvl1(time.Now().Format(time.RFC3339), "Computing SNP missingness filter")

	mpcPar := qc.general.mpcObj
	rtype := mpcPar[0].GetRType()
	pid := mpcPar[0].GetPid()
	gwasParams := qc.general.gwasParams
	numSnp := qc.general.gwasParams.numSnps
	useBoolean := mpcPar[0].GetBooleanShareFlag()

	totalInds := Sum(gwasParams.NumInds())

	xCount := make([]int, numSnp)

	// Take a pass over geno block files to get missingness information
	if pid > 0 {

		shifts := make([]int, len(qc.general.genoBlockSizes))
		for i := range shifts {
			if i > 0 {
				shifts[i] = shifts[i-1] + qc.general.genoBlockSizes[i-1]
			}
		}

		for i, genoFs := range qc.general.genoBlocks {

			for indiv, rowIndex := genoFs.NextRow(), 0; indiv != nil; indiv, rowIndex = genoFs.NextRow(), rowIndex+1 {
				for j, x := range indiv {
					snp := int(x)
					if snp >= 0 { // Not missing
						xCount[shifts[i]+j] += 1
					}
				}
			}

			genoFs.Reset()
		}

	}

	/* SNP missing rate */
	xCountRV := mpc_core.IntToRVec(rtype, xCount)
	lb := int((1 - qc.filterParams.GenoMissBound) * float64(totalInds))

	start := time.Now()
	log.LLvl1(time.Now().Format(time.RFC3339), "Secure comparison for SNP missingness ... ")
	gmissFilt := mpcPar.NotLessThanPublic(xCountRV, rtype.FromInt(lb), useBoolean) // Secure comparison
	log.LLvl1(time.Now().Format(time.RFC3339), "done. ", time.Since(start))

	snpFilt := mpcPar.RevealSymVec(gmissFilt)

	// Update geno file stream
	jkeep := make([]bool, len(snpFilt))
	if pid > 0 {
		for i := range jkeep {
			jkeep[i] = snpFilt[i].Uint64() != 0
		}
	}

	return jkeep
}

func (qc *QC) SNPMAFAndHWEFilters() []bool {

	log.LLvl1(time.Now().Format(time.RFC3339), "Computing SNP filters (minor allele frequency and Hardy-Weinberg equilibrium)")

	mpcPar := qc.general.mpcObj
	rtype := mpcPar[0].GetRType()
	pid := mpcPar[0].GetPid()
	fracBits := mpcPar[0].GetFracBits()
	numSnp := qc.filtNumSnps
	useBoolean := mpcPar[0].GetBooleanShareFlag()

	mafBound := qc.filterParams.MafLowerBound
	hweBound := rtype.FromFloat64(qc.filterParams.HweUpperBound, fracBits)

	xSum := make([]int, numSnp)
	xCount := make([]int, numSnp)

	// Over control cohort only
	xSumCtrl := make([]int, numSnp)
	xCountCtrl := make([]int, numSnp)
	genoObservedCtrl := make([][]int, 3)
	for i := range genoObservedCtrl {
		genoObservedCtrl[i] = make([]int, numSnp)
	}

	// Take a pass over geno block files to get missing and dosage information across dataset
	if pid > 0 {
		log.LLvl1(time.Now().Format(time.RFC3339), "Scanning the input files ... ")
		start := time.Now()

		phenoF := qc.general.pheno

		shifts := make([]int, len(qc.general.genoBlocks))
		for i := range shifts {
			if i > 0 {
				shifts[i] = shifts[i-1] + int(qc.general.genoBlocks[i-1].NumColsToKeep())
			}
		}

		for i, genoFs := range qc.general.genoBlocks {

			for indiv, rowIndex := genoFs.NextRow(), 0; indiv != nil; indiv, rowIndex = genoFs.NextRow(), rowIndex+1 {

				yi := int(phenoF.At(rowIndex, 0))

				for j, x := range indiv {
					snp := int(x)
					if snp >= 0 { // Not missing
						xSum[shifts[i]+j] += snp
						xCount[shifts[i]+j] += 2

						if yi < 1 { // Control cohort
							xSumCtrl[shifts[i]+j] += snp
							xCountCtrl[shifts[i]+j] += 2
							genoObservedCtrl[snp][shifts[i]+j]++
						}
					}
				}
			}

			genoFs.Reset()
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "done.", time.Since(start))
	}

	/* Minor allele frequency (MAF) */

	// Target: bound <= s/c <= 1 - bound
	// Compute as: (2*s-c)^2 <= c^2 * (2*bound - 1)^2
	xCountRV := mpc_core.IntToRVec(rtype, xCount)
	xSumRV := mpc_core.IntToRVec(rtype, xSum)

	{ // DEBUG
		if pid > 0 {
			SaveIntVectorToFile(qc.general.CachePath("gkeep_test_xcount.txt"), mpcPar.RevealSymVec(xCountRV).ToInt())
			SaveIntVectorToFile(qc.general.CachePath("gkeep_test_xsum.txt"), mpcPar.RevealSymVec(xSumRV).ToInt())
			log.LLvl1(time.Now().Format(time.RFC3339), "SNP XCount and XSum wrote to cache (for debugging)")
		}
	}

	if pid > 0 {
		xSumRV.MulScalar(rtype.FromInt(2))
		xSumRV.Sub(xCountRV)
	}
	xSumSq := mpcPar.SSMultElemVec(xSumRV, xSumRV)
	xCountSq := mpcPar.SSMultElemVec(xCountRV, xCountRV)

	if pid > 0 {
		prec := 20
		bound := rtype.FromFloat64(math.Pow(2*mafBound-1.0, 2), prec)

		xCountSq.MulScalar(bound)
		xSumSq.MulScalar(rtype.FromUint64(uint64(1) << prec))

		xCountSq.Sub(xSumSq)
	}

	start := time.Now()
	log.LLvl1(time.Now().Format(time.RFC3339), "IsPositive for MAF filter ... ")

	mafFilt := mpcPar.IsPositive(xCountSq, useBoolean)

	log.LLvl1(time.Now().Format(time.RFC3339), "done.", time.Since(start))

	{ // DEBUG
		if pid > 0 {
			SaveFloatVectorToFile(qc.general.CachePath("gkeep_test_ispos.txt"), mpcPar.RevealSymVec(xCountSq).ToFloat(20))

			mafFiltReveal := mpcPar.RevealSymVec(mafFilt)
			tmp := make([]bool, len(mafFiltReveal))
			for i := range tmp {
				tmp[i] = mafFiltReveal[i].Uint64() != 0
			}
			writeFilterToFile(qc.general.CachePath("gkeep_maf_only.txt"), tmp, false)
			log.LLvl1(time.Now().Format(time.RFC3339), "SNP miss filter wrote to cache:", qc.general.CachePath("gkeep_maf_only.txt"))
		}
	}

	/* Hardy-Weinberg equlibrium (over control cohort only) */

	xSumCtrlRV := mpc_core.IntToRVec(rtype, xSumCtrl)     // alpha
	xCountCtrlRV := mpc_core.IntToRVec(rtype, xCountCtrl) // beta
	genoObservedCtrlRV := mpc_core.IntToRMat(rtype, genoObservedCtrl)

	// Calculate expected genotype frequencies eAA, eAa, eaa
	// alpha = dosage sum over observed controls
	// beta = 2 * observed controls size
	// eaa = (beta - alpha)^2 / (2 * beta)
	// eAa = (2 * alpha * (beta - alpha)) / (2 * beta)
	// eAA = alpha^2 / (2 * beta)
	// Note: Denominator multiplied to the bound
	xCountCtrlConst := xCountCtrlRV.Copy()
	xCountCtrlConst.MulScalar(rtype.FromInt(2)) // 2 * beta
	xCountCtrlRV.Sub(xSumCtrlRV)                // beta - alpha

	expected := make([]mpc_core.RVec, 3)                       // AA, Aa, aa
	expected[2] = mpcPar.SSMultElemVec(xSumCtrlRV, xSumCtrlRV) // alpha^2
	xSumCtrlRV.MulScalar(rtype.FromInt(2))
	expected[1] = mpcPar.SSMultElemVec(xSumCtrlRV, xCountCtrlRV)   // 2 * alpha * (beta - alpha)
	expected[0] = mpcPar.SSMultElemVec(xCountCtrlRV, xCountCtrlRV) // (beta - alpha)^2

	// Compute chi-squared statistics
	// Sum { (genoObs * (2 * beta) - expected^2 / expected }
	chiSq := mpc_core.InitRVec(rtype.Zero(), numSnp)
	for i := range expected {
		tmp := mpcPar.SSMultElemVec(xCountCtrlConst, genoObservedCtrlRV[i])
		tmp.Sub(expected[i])
		tmp = mpcPar.SSMultElemVec(tmp, tmp)

		start = time.Now()

		divOut := mpcPar.Divide(tmp, expected[i], useBoolean)
		chiSq.Add(divOut)

		log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Division round %d/%d complete, %s", i+1, len(expected), time.Since(start)))
	}

	xCountCtrlConst.MulScalar(hweBound)

	start = time.Now()
	log.LLvl1(time.Now().Format(time.RFC3339), "Secure comparison for HWE ... ")
	hweFilt := mpcPar.LessThan(chiSq, xCountCtrlConst, useBoolean)
	log.LLvl1(time.Now().Format(time.RFC3339), "done.", time.Since(start))

	// Multiply all SNP filters and reveal
	snpFilt := mpcPar.SSMultElemVec(mafFilt, hweFilt)
	snpFilt = mpcPar.RevealSymVec(snpFilt)

	jkeep := make([]bool, len(snpFilt))
	if pid > 0 {
		for i := range jkeep {
			jkeep[i] = snpFilt[i].Uint64() != 0
		}
	}

	return jkeep
}

func (qc *QC) QualityControlProtocolWithPrecomputedGenoStats(useCache bool) {
	mpc := qc.general.mpcObj[0]
	pid := mpc.GetPid()

	var snpFilt []bool
	var nSnpFilt int
	snpFiltCache := qc.general.CachePath("gkeep.txt")
	if useCache {
		if pid > 0 {
			snpFilt = readFilterFromFile(snpFiltCache, qc.general.gwasParams.numSnps, false)
			log.LLvl1(time.Now().Format(time.RFC3339), "SNP QC filter loaded from cache:", snpFiltCache)
		}
	} else {

		genoStatsFile := qc.general.config.GenoCountFile
		var ac, gc [][]uint32
		var miss []uint32
		if pid > 0 {
			ac, gc, miss = ReadGenoStatsFromFile(genoStatsFile, qc.general.gwasParams.numSnps)
		} else { // pid == 0
			ac = make([][]uint32, 2)
			for i := range ac {
				ac[i] = make([]uint32, qc.general.gwasParams.numSnps)
			}
			gc = make([][]uint32, 3)
			for i := range gc {
				gc[i] = make([]uint32, qc.general.gwasParams.numSnps)
			}
			miss = make([]uint32, qc.general.gwasParams.numSnps)
		}

		snpFilt = qc.SNPFilterWithPrecomputedStats(ac, gc, miss, useCache)

		if pid > 0 {
			writeFilterToFile(snpFiltCache, snpFilt, false)
			log.LLvl1(time.Now().Format(time.RFC3339), "SNP QC filter wrote to cache:", snpFiltCache)
		}
	}

	if pid > 0 {
		nSnpFilt = SumBool(snpFilt)

		log.LLvl1(time.Now().Format(time.RFC3339), "Total # SNPs after QC filters (missingness/MAF/HWE):",
			qc.general.gwasParams.numSnps, "->", nSnpFilt)
	}

	// Share reduced SNP count with Party 0
	if pid == mpc.GetHubPid() {
		mpc.Network.SendInt(nSnpFilt, 0)
	} else if pid == 0 {
		nSnpFilt = mpc.Network.ReceiveInt(mpc.GetHubPid())
	}

	qc.filtNumSnps = nSnpFilt
	qc.filtNumInds = qc.general.gwasParams.numInds

	// Update gwasParams
	qc.general.gwasParams.SetFiltCounts(qc.filtNumInds, qc.filtNumSnps)
	qc.general.gwasParams.SetSnpFilt(snpFilt)

	return
}

//QualityControlProtocol (1) applies SNP filters and individual filters to input data, (2) filters input wrt filters
func (qc *QC) QualityControlProtocol(useCache bool) {
	mpc := qc.general.mpcObj[0]
	pid := mpc.GetPid()

	var snpMissFilt []bool
	var nSnpMissFilt int
	snpMissCache := qc.general.CachePath("gkeep_miss.txt")
	if useCache {
		if pid > 0 {
			snpMissFilt = readFilterFromFile(snpMissCache, qc.general.gwasParams.numSnps, false)
			log.LLvl1(time.Now().Format(time.RFC3339), "SNP miss filter loaded from cache:", snpMissCache)
		}
	} else {

		snpMissFilt = qc.SNPMissFilter()

		if pid > 0 {
			writeFilterToFile(snpMissCache, snpMissFilt, false)
			log.LLvl1(time.Now().Format(time.RFC3339), "SNP miss filter wrote to cache:", snpMissCache)
		}
	}

	if pid > 0 {
		// Update geno streams
		shift := 0
		for _, geno := range qc.general.genoBlocks {
			m := int(geno.NumCols())
			geno.UpdateColFilt(snpMissFilt[shift : shift+m])
			shift += m
		}
		nSnpMissFilt = SumBool(snpMissFilt)

		log.LLvl1(time.Now().Format(time.RFC3339), "# SNPs after missingness filter:",
			qc.general.gwasParams.numSnps, "->", nSnpMissFilt)

		qc.filtNumSnps = nSnpMissFilt
	}

	var indFilt []bool
	var nIndFilt int
	indFiltCache := qc.general.CachePath("ikeep.txt")
	if pid > 0 {
		if useCache {
			indFilt = readFilterFromFile(indFiltCache, qc.general.gwasParams.numInds[pid], false)
			log.LLvl1(time.Now().Format(time.RFC3339), "Individual filter loaded from cache:", indFiltCache)
		} else {
			indFilt = qc.IndividualMissAndHetFilters()

			writeFilterToFile(indFiltCache, indFilt, false)
			log.LLvl1(time.Now().Format(time.RFC3339), "Individual filter wrote to cache:", indFiltCache)
		}

		// Update geno streams
		for _, genoFs := range qc.general.genoBlocks {
			genoFs.UpdateRowFilt(indFilt)
		}
		nIndFilt = SumBool(indFilt)

		log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Number of individuals: %d -> %d",
			qc.general.gwasParams.numInds[pid], nIndFilt))
	}

	// Share filtered snp/individual count with other parties
	var vec []uint64
	if pid == mpc.GetHubPid() {
		vec = make([]uint64, mpc.GetNParty())
		vec[pid] = uint64(nIndFilt)
		for other := 1; other < mpc.GetNParty(); other++ {
			if other != pid {
				vec[other] = uint64(mpc.Network.ReceiveInt(other))
			}
		}
		for other := 0; other < mpc.GetNParty(); other++ { // Include Party 0!
			if other != pid {
				mpc.Network.SendIntVector(vec, other)
			}
		}

		mpc.Network.SendInt(qc.filtNumSnps, 0)
	} else if pid > 0 {
		mpc.Network.SendInt(nIndFilt, mpc.GetHubPid())
		vec = mpc.Network.ReceiveIntVector(mpc.GetNParty(), mpc.GetHubPid())
	} else if pid == 0 {
		vec = mpc.Network.ReceiveIntVector(mpc.GetNParty(), mpc.GetHubPid())
		nSnpMissFilt = mpc.Network.ReceiveInt(mpc.GetHubPid())
		qc.filtNumSnps = nSnpMissFilt
	}

	for i := range vec {
		qc.filtNumInds[i] = int(vec[i])
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "# SNPs before maf/hwe filters:", qc.filtNumSnps)
	log.LLvl1(time.Now().Format(time.RFC3339), "# Individuals:", qc.filtNumInds[1:])

	// Compute SNP filters
	var snpFilt []bool
	var nSnpFilt int
	snpMafHweCache := qc.general.CachePath("gkeep_maf_hwe.txt")
	if useCache {
		if pid > 0 {
			snpFilt = readFilterFromFile(snpMafHweCache, qc.filtNumSnps, false)
			log.LLvl1(time.Now().Format(time.RFC3339), "SNP MAF/HWE filter loaded from cache:", snpMafHweCache)
		}
	} else {
		snpFilt = qc.SNPMAFAndHWEFilters()

		if pid > 0 {
			writeFilterToFile(snpMafHweCache, snpFilt, false)
			log.LLvl1(time.Now().Format(time.RFC3339), "SNP MAF/HWE filter wrote to cache:", snpMafHweCache)
		}
	}

	// Update geno streams
	if pid > 0 {
		shift := 0
		for _, geno := range qc.general.genoBlocks {
			m := int(geno.NumColsToKeep())
			geno.UpdateColFilt(snpFilt[shift : shift+m])
			shift += m
		}
		nSnpFilt = SumBool(snpFilt)
	}

	// Share reduced SNP count with Party 0
	if pid == mpc.GetHubPid() {
		mpc.Network.SendInt(nSnpFilt, 0)
	} else if pid == 0 {
		nSnpFilt = mpc.Network.ReceiveInt(mpc.GetHubPid())
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "# SNPs after MAF/HWE filters:",
		nSnpMissFilt, "->", nSnpFilt)

	qc.filtNumSnps = nSnpFilt

	// Update gwasParams
	qc.general.gwasParams.SetFiltCounts(qc.filtNumInds, qc.filtNumSnps)
	qc.general.gwasParams.SetSnpFilt(snpFilt)

	// Filter pheno and cov data
	if pid > 0 {
		qc.general.pheno = FilterVec(qc.general.pheno, indFilt)
		qc.general.cov = FilterMat(qc.general.cov, OnesBool(qc.general.gwasParams.NumCov()), indFilt)
	}

	return
}

func FilterMat(X *mat.Dense, colkeep, rowkeep []bool) *mat.Dense {
	r, c := X.Dims()
	xdata := X.RawMatrix().Data

	rnew, cnew := 0, 0
	for i := range rowkeep {
		if rowkeep[i] {
			rnew++
		}
	}
	for i := range colkeep {
		if colkeep[i] {
			cnew++
		}
	}

	xdataNew := make([]float64, rnew*cnew)
	index := 0

	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			if rowkeep[i] && colkeep[j] {
				xdataNew[index] = xdata[(c*i)+j]
				index++
			}

		}
	}

	Xnew := mat.NewDense(rnew, cnew, xdataNew)

	return Xnew
}

func FilterVec(vec *mat.Dense, keep []bool) *mat.Dense {
	n, _ := vec.Dims()
	vecData := vec.RawMatrix().Data

	newlen := 0
	for i := range keep {
		if keep[i] {
			newlen++
		}
	}
	vecDataNew := make([]float64, newlen)
	index := 0
	for i := 0; i < n; i++ {
		if keep[i] {
			vecDataNew[index] = vecData[i]
			index++
		}

	}

	vecNew := mat.NewDense(newlen, 1, vecDataNew)

	return vecNew
}
