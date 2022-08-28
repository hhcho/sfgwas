package gwas

type FilterParams struct {
	MafLowerBound float64
	HweUpperBound float64
	GenoMissBound float64
	IndMissBound  float64
	HetLowerBound float64
	HetUpperBound float64
}

type GWASParams struct {
	numInds    []int
	numCovs    int
	numSnps    int
	numSnpsPCA int
	numPCs     int

	numFiltInds []int
	numFiltSnps int

	snpFilt []bool

	minSnpDist uint64

	skipQC bool
	runPCA bool // LMM or PCA
}

func InitGWASParams(numInds []int, numSnps, numCovs, numPCs, minSnpDist int) *GWASParams {
	gwasParams := &GWASParams{
		numInds:    numInds,
		numSnps:    numSnps,
		numCovs:    numCovs,
		numPCs:     numPCs,
		minSnpDist: uint64(minSnpDist),
	}
	return gwasParams
}

func (gwasParams *GWASParams) SnpFilt() []bool {
	return gwasParams.snpFilt
}

func (gwasParams *GWASParams) SetSnpFilt(filt []bool) {
	gwasParams.snpFilt = filt
}

func (gwasParams *GWASParams) SetPopStratMethod(s bool) {
	gwasParams.runPCA = s
}

func (gwasParams *GWASParams) GetPopStratMethod() bool {
	return gwasParams.runPCA
}

func (gwasParams *GWASParams) NumInds() []int {
	return gwasParams.numInds
}

func (gwasParams *GWASParams) NumSNP() int {
	return gwasParams.numSnps
}

func (gwasParams *GWASParams) NumCov() int {
	return gwasParams.numCovs
}

func (gwasParams *GWASParams) NumPC() int {
	return gwasParams.numPCs
}

func (gwasParams *GWASParams) SetNumPC(numPCs int) {
	gwasParams.numPCs = numPCs
}

func (gwasParams *GWASParams) MinSnpDistThreshold() uint64 {
	return gwasParams.minSnpDist
}

func (gwasParams *GWASParams) SetFiltCounts(filtInds []int, filtSnps int) {
	gwasParams.numFiltInds = filtInds
	gwasParams.numFiltSnps = filtSnps
}

func (gwasParams *GWASParams) SetNumSnpsPCA(numSnps int) {
	gwasParams.numSnpsPCA = numSnps
}

func (gwasParams *GWASParams) FiltNumInds() []int {
	return gwasParams.numFiltInds
}

func (gwasParams *GWASParams) FiltNumSNP() int {
	return gwasParams.numFiltSnps
}

func InitFilteringSettings(maflb, hwe, gmiss, imiss, hetlb, hetub float64) *FilterParams {
	filterParams := &FilterParams{
		MafLowerBound: maflb,
		HweUpperBound: hwe,
		GenoMissBound: gmiss,
		IndMissBound:  imiss,
		HetLowerBound: hetlb,
		HetUpperBound: hetub,
	}
	return filterParams
}
