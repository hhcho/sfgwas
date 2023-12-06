package mpc

import (
	"encoding/gob"
	"fmt"
	"math"
	"math/big"
	"time"

	"go.dedis.ch/onet/v3/log"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/ldsec/lattigo/v2/dckks"
	"github.com/ldsec/lattigo/v2/ring"

	"sync"
)

type MPC struct {
	rtype            mpc_core.RElem
	dataBits         int
	fracBits         int
	useBooleanShares bool
	Network          *Network

	divSqrtMaxLen int

	orLagrangeCache map[TypedKey]mpc_core.RVec
	lagrangeCache   map[TypedKey]mpc_core.RMat
	invPowCache     map[TypedKey]mpc_core.RElem
	pascalCache     map[TypedKey]mpc_core.RMat

	syncCounter int
}

type ParallelMPC []*MPC // Holds mpc environment for each thread for parallelization
type ParallelNetworks []*Network

func InitParallelMPCEnv(netObjs []*Network, rtype mpc_core.RElem, dataBits, fracBits int) []*MPC {

	mpcEnvParallel := make([]*MPC, len(netObjs))
	mpcEnvParallel[0] = initMPCEnv(netObjs[0], rtype, dataBits, fracBits)
	lagrangeCache := mpcEnvParallel[0].lagrangeCache

	for i := 1; i < len(netObjs); i++ {
		mpcEnvParallel[i] = &MPC{
			Network:          netObjs[i],
			dataBits:         dataBits,
			fracBits:         fracBits,
			useBooleanShares: true,
			orLagrangeCache:  make(map[TypedKey]mpc_core.RVec),
			lagrangeCache:    lagrangeCache, // reuse lagrange cache for all mpc instances
			invPowCache:      make(map[TypedKey]mpc_core.RElem),
			pascalCache:      make(map[TypedKey]mpc_core.RMat),
			syncCounter:      0,
			rtype:            rtype,
		}
	}

	return mpcEnvParallel
}

func initMPCEnv(netObj *Network, rtype mpc_core.RElem, dataBits, fracBits int) *MPC {
	gob.Register(mpc_core.LElem256Zero)
	gob.Register(mpc_core.LElem128Zero)
	gob.Register(mpc_core.LElem2N(0))
	gob.Register(mpc_core.LElemP(0))
	gob.Register(mpc_core.SElemDS(0))
	gob.Register(mpc_core.SElemC(0))
	gob.Register(mpc_core.RVec{nil})
	gob.Register(mpc_core.RMat{nil})

	mpcObj := &MPC{
		Network:         netObj,
		dataBits:        dataBits,
		fracBits:        fracBits,
		orLagrangeCache: make(map[TypedKey]mpc_core.RVec),
		lagrangeCache:   make(map[TypedKey]mpc_core.RMat),
		invPowCache:     make(map[TypedKey]mpc_core.RElem),
		pascalCache:     make(map[TypedKey]mpc_core.RMat),
		syncCounter:     0,
		rtype:           rtype,
	}

	mpcObj.InitLagrangeCache()
	return mpcObj
}

type TypedKey struct {
	key    int
	typeId uint8
}

func pascalmat(rtype mpc_core.RElem, pow int) mpc_core.RMat {
	res := make(mpc_core.RMat, pow+1)

	for i := range res {
		res[i] = make(mpc_core.RVec, pow+1)
		for j := range res[i] {
			if j > i {
				res[i][j] = rtype.Zero()
			} else if j == 0 || j == i {
				res[i][j] = rtype.One()
			} else {
				res[i][j] = res[i-1][j-1].Add(res[i-1][j])
			}
		}
	}
	return res
}

//GetPascalMat retrieves the pascal matrix associated with the given power
func (mpcObj *MPC) GetPascalMatrix(rtype mpc_core.RElem, pow int) mpc_core.RMat {
	k := TypedKey{pow, rtype.TypeID()}

	res, found := mpcObj.pascalCache[k]

	if found {
		// log.LLvl1(time.Now().Format(time.RFC3339), "status of cache: ", pascalCache)
		return res
	}

	//compute, update and return
	res = pascalmat(rtype, pow)
	mpcObj.pascalCache[k] = res
	return res
}

func (mpcObj *MPC) AssertSync() {
	pid := mpcObj.GetPid()
	check := mpcObj.syncCounter

	if pid == mpcObj.GetHubPid() {
		for other := 0; other < mpcObj.GetNParty(); other++ {
			if other == pid {
				continue
			}

			if check != mpcObj.Network.ReceiveInt(other) {
				log.Fatal("AssertSync counter check failed between parties", pid, "and", other)
			}
		}
	} else {
		mpcObj.Network.SendInt(check, mpcObj.GetHubPid())
	}

	for other := 0; other < mpcObj.GetNParty(); other++ {
		if other == pid {
			continue
		}

		mpcObj.Network.Rand.SwitchPRG(other)
		rCheck := int(mpcObj.Network.Rand.RandElem(mpcObj.GetRType()).Uint64())
		mpcObj.Network.Rand.RestorePRG()

		var otherCheck int

		if pid < other {
			mpcObj.Network.SendInt(rCheck, other)
			otherCheck = mpcObj.Network.ReceiveInt(other)
		} else {
			otherCheck = mpcObj.Network.ReceiveInt(other)
			mpcObj.Network.SendInt(rCheck, other)
		}

		if rCheck != otherCheck {
			log.Fatal("AssertSync PRG check failed between parties", pid, "and", other, ":", rCheck, "!=", otherCheck)
		}

	}

	mpcObj.syncCounter++
}

func (mpcObj *MPC) SetPid(p int) {
	mpcObj.Network.pid = p
}

func (mpcObj *MPC) GetPid() int {
	return mpcObj.Network.pid
}

func (mpcObj *MPC) SetHubPid(p int) {
	mpcObj.Network.hubPid = p
}

func (mpcObj *MPC) GetBooleanShareFlag() bool {
	return mpcObj.useBooleanShares
}

func (mpcObj *MPC) SetBooleanShareFlag(flag bool) {
	mpcObj.useBooleanShares = flag
}

func (mpcObj *MPC) SetDivSqrtMaxLen(num int) {
	mpcObj.divSqrtMaxLen = num
}

func (mpcObj *MPC) GetHubPid() int {
	return mpcObj.Network.hubPid
}

func (mpcObj *MPC) SetNParty(np int) {
	mpcObj.Network.NumParties = np
}

func (mpcObj *MPC) GetNParty() int {
	return mpcObj.Network.NumParties
}

func (mpcObj *MPC) GetDataBits() int {
	return mpcObj.dataBits
}

func (mpcObj *MPC) GetFracBits() int {
	return mpcObj.fracBits
}

func (mpcObj *MPC) SetFracBits(f int) {
	mpcObj.fracBits = f
}

func (mpcObj *MPC) GetRType() mpc_core.RElem {
	return mpcObj.rtype
}

func (mpcObj *MPC) GetMHEContext() *dckks.Context {
	return mpcObj.Network.dckksContext
}

func (mpcObj *MPC) GetCRPGen() *ring.UniformSampler {
	return mpcObj.Network.crpGen
}

func (mpcObj *MPC) InitLagrangeCache() {
	lagrangeCache := make(map[TypedKey]mpc_core.RMat)

	tableList := make([]mpc_core.RMat, 2)
	rtypeInputList := make([]mpc_core.RElem, 2)
	var rtype mpc_core.RElem
	var rtypeInput mpc_core.RElem
	var table mpc_core.RMat

	// Table 0: mpcObj.IsPositive
	rtype = mpcObj.rtype

	rtypeInput = mpc_core.SElemC(0)
	table = mpc_core.InitRMat(rtype.Zero(), 1, 2)
	if mpcObj.Network.pid > 0 {
		table[0][0] = rtype.One()
		table[0][1] = rtype.Zero()
	}
	tableList[0] = table
	rtypeInputList[0] = rtypeInput

	// Table 1: NormalizerEvenExp
	rtype = mpcObj.rtype

	rtypeInput = mpc_core.SElemDS(0)
	halfLen := mpcObj.dataBits / 2 // TODO: make this a constant
	table = mpc_core.InitRMat(rtype.Zero(), 2, halfLen+1)
	if mpcObj.Network.pid > 0 {
		for i := 0; i < halfLen+1; i++ {
			if i == 0 {
				table[0][i] = rtype.One()
				table[1][i] = rtype.One()
			} else {
				table[0][i] = table[0][i-1].Mul(rtype.FromInt(2))
				table[1][i] = table[1][i-1].Mul(rtype.FromInt(4))
			}
		}
	}
	tableList[1] = table
	rtypeInputList[1] = rtypeInput

	for t := range tableList {
		rtype := tableList[t].Type()
		rtypeInput := rtypeInputList[t]
		nrow, ncol := tableList[t].Dims()
		coeff := mpc_core.InitRMat(rtype.Zero(), nrow, (mpcObj.Network.NumParties-1)*ncol)
		modulus := rtype.FromUint64(rtypeInput.Modulus().Uint64())
		if mpcObj.Network.pid > 0 {
			for i := 0; i < nrow; i++ {
				x := make(mpc_core.RVec, (mpcObj.Network.NumParties-1)*ncol)
				y := make(mpc_core.RVec, (mpcObj.Network.NumParties-1)*ncol)
				for j := 0; j < ncol; j++ {
					val := tableList[t][i][j]
					x[j] = rtype.FromInt(j + 1)
					y[j] = val
					for p := 1; p < mpcObj.Network.NumParties-1; p++ {
						x[j+ncol*p] = x[j+ncol*(p-1)].Add(modulus)
						y[j+ncol*p] = val
					}
				}

				//fmt.Println("Generating table", t, "row", i)
				//fmt.Println("x:", x)
				//fmt.Println("y:", y)

				coeff[i] = lagrangeInterp(x, y)
			}
		}

		key := TypedKey{t, rtype.TypeID()}
		lagrangeCache[key] = coeff
	}

	mpcObj.lagrangeCache = lagrangeCache
}

// Interpolate (x_1,y_1), (x_2,y_2), ... (x_n,y_n)
// Return polynomial coefficients in the order of (1, x, x^2, ..., x^(n-1))
func lagrangeInterp(x, y mpc_core.RVec) mpc_core.RVec {
	if len(x) != len(y) {
		panic("Length of vectors are mismatching")
	}

	n := len(y)
	rtype := x[0]

	invTable := make(map[uint64]mpc_core.RElem)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != j {
				diff := x[j].Sub(x[i])

				_, flag := invTable[diff.Uint64()]
				if !flag {
					invTable[diff.Uint64()] = diff.Inv()
				}
			}
		}
	}

	//initalize denom and numer
	numer := mpc_core.InitRMat(rtype.Zero(), n, n)
	numer[0] = y.Copy()
	denomInv := mpc_core.InitRVec(rtype.One(), n)

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != j {

				for k := n - 1; k >= 0; k-- {
					if k == 0 {
						numer[k][j] = rtype.Zero().Sub(x[i].Mul(numer[k][j]))
					} else {
						numer[k][j] = numer[k-1][j].Sub(x[i].Mul(numer[k][j]))
					}
				}

				denomInv[i] = denomInv[i].Mul(invTable[x[j].Sub(x[i]).Uint64()].Neg())
			}
		}
	}

	out := mpc_core.RMultMatVec(numer, denomInv)

	return out
}

// Interpolate (1,y_1), (2,y_2), ... (n,y_n)
// Return polynomial coefficients in the order of (1, x, x^2, ..., x^(n-1))
func lagrangeInterpSimple(y mpc_core.RVec) mpc_core.RVec {
	rtype := y[0]
	n := len(y)
	x := make(mpc_core.RVec, n)
	for i := range x {
		x[i] = rtype.FromInt(i + 1)
	}
	return lagrangeInterp(x, y)
}

func (mpcObj *MPC) RevealSym(a mpc_core.RElem) mpc_core.RElem {
	pid := mpcObj.Network.pid
	if pid == 0 {
		return a
	}
	ar := mpcObj.RevealSymMat(mpc_core.RMat{mpc_core.RVec{a}})
	return ar[0][0]
}

func (mpcObj *MPC) RevealSymVec(a mpc_core.RVec) mpc_core.RVec {
	pid := mpcObj.Network.pid
	if pid == 0 {
		return a
	}
	ar := mpcObj.RevealSymMat(mpc_core.RMat{a})
	return ar[0]
}

func (mpcObj *MPC) RevealSymMat(a mpc_core.RMat) mpc_core.RMat {
	pid := mpcObj.Network.pid
	if pid == 0 {
		return a
	}

	rtype := a.Type()
	nr, nc := a.Dims()

	ar := a.Copy()

	for p := 1; p < mpcObj.Network.NumParties; p++ {
		if p < pid {
			mpcObj.Network.SendRData(a, p)
			r := mpcObj.Network.ReceiveRMat(rtype, nr, nc, p)
			ar.Add(r)

			//DEBUG
			//fmt.Println(pid, "send", p, a[0][:5])
			//fmt.Println(pid, "receive", p, r[0][:5])
		} else if p > pid {
			r := mpcObj.Network.ReceiveRMat(rtype, nr, nc, p)
			ar.Add(r)
			mpcObj.Network.SendRData(a, p)

			//DEBUG
			//fmt.Println(pid, "send", p, a[0][:5])
			//fmt.Println(pid, "receive", p, r[0][:5])
		}
	}

	return ar
}

// The i-th column of output contains 0,1,...,pow powers of a[i]
func (mpcObj *MPC) Powers(a mpc_core.RVec, pow int) mpc_core.RMat {
	pid := mpcObj.Network.pid
	n := len(a)
	rtype := a.Type()

	if pow == 1 {
		b := mpc_core.InitRMat(rtype.Zero(), 2, n)
		if pid > 0 {
			if pid == 1 {
				b[0] = mpc_core.InitRVec(rtype.One(), n)
			}
			b[1] = a.Copy()
		}
		return b
	}

	// Pow > 1
	ar, am := mpcObj.BeaverPartitionVec(a)

	if pid == 0 {
		ampow := make(mpc_core.RMat, pow-1)
		ampow[0] = am.Copy()
		ampow[0].MulElem(am)
		for p := 1; p < pow-1; p++ {
			ampow[p] = ampow[p-1].Copy()
			ampow[p].MulElem(am)
		}

		// Secret share ampow
		for p := 1; p < mpcObj.Network.NumParties-1; p++ {
			mpcObj.Network.Rand.SwitchPRG(p)
			share := mpcObj.Network.Rand.RandMat(rtype, pow-1, n)
			mpcObj.Network.Rand.RestorePRG()
			ampow.Sub(share)
		}

		mpcObj.Network.SendRData(ampow, mpcObj.Network.NumParties-1)

		b := mpc_core.InitRMat(rtype.Zero(), pow+1, n)
		return b
	}
	var ampow mpc_core.RMat
	if pid == mpcObj.Network.NumParties-1 {
		ampow = mpcObj.Network.ReceiveRMat(rtype, pow-1, n, 0)
	} else {
		mpcObj.Network.Rand.SwitchPRG(0)
		ampow = mpcObj.Network.Rand.RandMat(rtype, pow-1, n)
		mpcObj.Network.Rand.RestorePRG()
	}

	// Compute powers of ar
	arpow := make(mpc_core.RMat, pow-1)
	arpow[0] = ar.Copy()
	arpow[0].MulElem(ar)
	for p := 1; p < pow-1; p++ {
		arpow[p] = arpow[p-1].Copy()
		arpow[p].MulElem(ar)
	}

	var t mpc_core.RMat
	t = mpcObj.GetPascalMatrix(rtype, pow)

	b := mpc_core.InitRMat(rtype.Zero(), pow+1, n)
	if pid == 1 {
		b[0] = mpc_core.InitRVec(rtype.One(), n)
	}
	b[1] = a.Copy()

	for p := 2; p <= pow; p++ {
		if pid == 1 {
			b[p] = arpow[p-2].Copy()
		}

		if p == 2 {
			for i := range b[p] {
				b[p][i] = b[p][i].Add(t[p][1].Mul(ar[i].Mul(am[i])))
			}
		} else {
			for i := range b[p] {
				b[p][i] = b[p][i].Add(t[p][1].Mul(arpow[p-3][i].Mul(am[i])))
			}

			for j := 2; j <= p-2; j++ {
				for i := range b[p] {
					b[p][i] = b[p][i].Add(t[p][j].Mul(arpow[p-2-j][i].Mul(ampow[j-2][i])))
				}
			}

			for i := range b[p] {
				b[p][i] = b[p][i].Add(t[p][p-1].Mul(ar[i].Mul(ampow[p-3][i])))
			}
		}

		b[p].Add(ampow[p-2])
	}

	return b
}

func (mpcObj *MPC) EvaluatePoly(a mpc_core.RVec, coeff mpc_core.RMat) mpc_core.RMat {
	pid := mpcObj.Network.pid
	n := len(a)
	rtype := a.Type()
	npoly, deg := coeff.Dims()
	deg -= 1

	apow := mpcObj.Powers(a, deg)

	if pid == 0 {
		return mpc_core.InitRMat(rtype.Zero(), npoly, n)
	}
	return mpc_core.RMultMat(coeff, apow)
}

func (mpcObj *MPC) FanInOr(a mpc_core.RMat) mpc_core.RVec {
	pid := mpcObj.Network.pid
	_, d := a.Dims()
	rtype := a.Type()

	rowSum := a.Sum(0)
	if pid == 1 {
		rowSum.AddScalar(rtype.One())
	}

	key := TypedKey{d + 1, rtype.TypeID()}
	coeff, found := mpcObj.orLagrangeCache[key]
	if !found {
		y := mpc_core.InitRVec(rtype.One(), d+1)
		y[0] = rtype.Zero()
		coeff = lagrangeInterpSimple(y)
		mpcObj.orLagrangeCache[key] = coeff
	}

	return mpcObj.EvaluatePoly(rowSum, mpc_core.RMat{coeff})[0]
}

// Given a vector of numbers with bit length nbits,
// extract bits from each number and output as a matrix of rtypeBit elements
// If rtypeBit is BElem, then bit packing is performed
func numToBits(rtypeBit mpc_core.RElem, a mpc_core.RVec, nbits int) mpc_core.RMat {
	n := len(a)

	isBElem := rtypeBit.TypeID() == mpc_core.BElemUniqueID
	nelem := 1 + ((nbits - 1) / mpc_core.BElemBits)

	out := make(mpc_core.RMat, n)
	for i := range out {
		if isBElem {
			out[i] = mpc_core.InitRVec(rtypeBit.Zero(), nelem)
		} else {
			out[i] = make(mpc_core.RVec, nbits)
		}

		for j := 0; j < nbits; j++ {
			bit := a[i].GetBit(nbits - 1 - j)
			if isBElem {
				index := nelem - 1 - ((nbits - 1 - j) / mpc_core.BElemBits)
				bitPos := (nbits - 1 - j) % mpc_core.BElemBits
				out[i][index] = out[i][index].Add(mpc_core.BElem(uint64(bit) << bitPos))
			} else {
				if bit > 0 {
					out[i][j] = rtypeBit.One()
				} else {
					out[i][j] = rtypeBit.Zero()
				}
			}
		}
	}
	return out
}

// Secret share random numbers in rtype (output as RVec)
// and bit-wise secret shares of the same numbers in rtypeBit (output as RMat)
// Party 0 obtains underlying values
// If bitSample is true, uniformly samples from [0, 2^nbitsSample);
// otherwise samples from [0, rtype.Modulus())
func (mpcObj *MPC) ShareRandomBits(rtype, rtypeBit mpc_core.RElem, n, nbitsOut, nbitsSample int, bitSample bool) (mpc_core.RVec, mpc_core.RMat) {
	pid := mpcObj.Network.pid

	isBElem := rtypeBit.TypeID() == mpc_core.BElemUniqueID
	nelem := nbitsOut
	if isBElem {
		nelem = 1 + (nbitsOut-1)/mpc_core.BElemBits
	}

	if pid == 0 {
		var r mpc_core.RVec
		if bitSample {
			r = mpcObj.Network.Rand.RandVecBits(rtype, n, nbitsSample)
		} else {
			r = mpcObj.Network.Rand.RandVec(rtype, n)
		}
		rBits := numToBits(rtypeBit, r, nbitsOut)

		// Secret share r and rBits
		for p := 1; p < mpcObj.Network.NumParties-1; p++ {
			mpcObj.Network.Rand.SwitchPRG(p)
			mask := mpcObj.Network.Rand.RandVec(rtype, n)
			maskBits := mpcObj.Network.Rand.RandMat(rtypeBit, n, nelem)
			mpcObj.Network.Rand.RestorePRG()

			r.Sub(mask)
			rBits.Sub(maskBits)
		}

		mpcObj.Network.SendRData(r, mpcObj.Network.NumParties-1)
		mpcObj.Network.SendRData(rBits, mpcObj.Network.NumParties-1)

		return r, rBits
	}

	var r mpc_core.RVec
	var rBits mpc_core.RMat
	if pid == mpcObj.Network.NumParties-1 {
		r = mpcObj.Network.ReceiveRVec(rtype, n, 0)
		rBits = mpcObj.Network.ReceiveRMat(rtypeBit, n, nelem, 0)
	} else {
		mpcObj.Network.Rand.SwitchPRG(0)
		r = mpcObj.Network.Rand.RandVec(rtype, n)
		rBits = mpcObj.Network.Rand.RandMat(rtypeBit, n, nelem)
		mpcObj.Network.Rand.RestorePRG()
	}

	return r, rBits
}

// Return prefix-OR of each row in a containing secret shares of bits
func (mpcObj *MPC) PrefixOr(a mpc_core.RMat) mpc_core.RMat {
	pid := mpcObj.Network.pid
	n, k := a.Dims()
	rtype := a.Type()

	// Find the next largest square integer
	L := int(math.Ceil(math.Sqrt(float64(k))))
	L2 := L * L

	// TODO: add sanity check
	//if a.Modulus() <= L + 1 {
	//	panic("Modulus too small for PrefixOr")
	//}

	// Zero pad to L2 bits, then reshape to n * L by L matrix
	aPadded := mpc_core.InitRMat(rtype.Zero(), n*L, L)
	if pid > 0 {
		for i := 0; i < n; i++ {
			for s1 := 0; s1 < L; s1++ {
				for s2 := 0; s2 < L; s2++ {
					j := s1*L + s2
					if j < L2-k {
						aPadded[L*i+s1][s2] = rtype.Zero()
					} else {
						aPadded[L*i+s1][s2] = a[i][j-L2+k]
					}
				}
			}
		}
	}

	x := mpcObj.FanInOr(aPadded)

	xPre := mpc_core.InitRMat(rtype.Zero(), n*L, L)
	if pid > 0 {
		for i := 0; i < n; i++ {
			for s1 := 0; s1 < L; s1++ {
				xpi := L*i + s1
				for s2 := 0; s2 < L; s2++ {
					if s2 <= s1 {
						xPre[xpi][s2] = x[L*i+s2]
					} // 0 otherwise
				}
			}
		}
	}

	y := mpcObj.FanInOr(xPre)
	xPre = nil

	f := mpc_core.InitRMat(rtype.Zero(), n, L)

	if pid > 0 {
		for i := 0; i < n; i++ {
			for j := 0; j < L; j++ {
				if j == 0 {
					f[i][j] = x[L*i]
				} else {
					f[i][j] = y[L*i+j].Sub(y[L*i+j-1])
				}
			}
		}
	}
	x = nil

	fr, fm := mpcObj.BeaverPartitionMat(f)
	apr, apm := mpcObj.BeaverPartitionMat(aPadded)

	rmat := make(mpc_core.RMat, L)
	mmat := make(mpc_core.RMat, L)
	c := make(mpc_core.RMat, n)
	for i := range c {
		// 1xL times LxL matrix multiplication
		for j := 0; j < L; j++ {
			rmat[j] = apr[L*i+j]
			mmat[j] = apm[L*i+j]
		}
		c[i] = mpcObj.BeaverMultMat(mpc_core.RMat{fr[i]}, mpc_core.RMat{fm[i]}, rmat, mmat)[0]
	}

	c = mpcObj.BeaverReconstructMat(c)
	aPadded, apr, apm = nil, nil, nil
	rmat, mmat = nil, nil

	cPre := mpc_core.InitRMat(rtype.Zero(), n*L, L)
	if pid > 0 {
		for i := 0; i < n; i++ {
			for s1 := 0; s1 < L; s1++ {
				cpi := L*i + s1
				for s2 := 0; s2 < L; s2++ {
					if s2 <= s1 {
						cPre[cpi][s2] = c[i][s2]
					} // 0 otherwise
				}
			}
		}
	}
	c = nil

	bdotVec := mpcObj.FanInOr(cPre)
	cPre = nil

	bdot := mpc_core.InitRMat(rtype.Zero(), n, L)
	if pid > 0 {
		for i := 0; i < n; i++ {
			for j := 0; j < L; j++ {
				bdot[i][j] = bdotVec[L*i+j]
			}
		}
	}
	bdotVec = nil

	bdotr, bdotm := mpcObj.BeaverPartitionMat(bdot)
	s := make(mpc_core.RMat, n*L) // Concatenation of n LxL matrices
	for i := 0; i < n; i++ {
		// Take two length-L vectors, output LxL outer product matrix
		tmp := mpcObj.BeaverMultMat(mpc_core.RMat{fr[i]}.Transpose(), mpc_core.RMat{fm[i]}.Transpose(), mpc_core.RMat{bdotr[i]}, mpc_core.RMat{bdotm[i]})
		for j := 0; j < L; j++ {
			s[L*i+j] = tmp[j]
		}
	}
	s = mpcObj.BeaverReconstructMat(s)
	bdot, bdotr, bdotm = nil, nil, nil

	b := mpc_core.InitRMat(rtype.Zero(), n, k)
	if pid > 0 {
		for i := 0; i < n; i++ {
			for j := 0; j < k; j++ {
				jPad := L2 - k + j
				il := int(jPad / L)
				jl := jPad - il*L
				b[i][j] = s[L*i+il][jl].Add(y[L*i+il]).Sub(f[i][il])
			}
		}
	}

	return b
}

func (mpcObj *MPC) TableLookup(a mpc_core.RVec, tableId int) mpc_core.RMat {
	rtype := a.Type()

	key := TypedKey{tableId, rtype.TypeID()}
	coeff, found := mpcObj.lagrangeCache[key]
	if !found {
		panic(fmt.Sprintf("Requested table for TableLookup is invalid: tableId %d typeID %d", key.key, key.typeId))
	}

	return mpcObj.EvaluatePoly(a, coeff)
}

func (mpcObj *MPC) TableLookupWithShareConversion(rtypeOut mpc_core.RElem, a mpc_core.RVec, tableId int) mpc_core.RMat {
	n := len(a)
	aConv := make(mpc_core.RVec, n)
	for i := range aConv {
		aConv[i] = rtypeOut.FromUint64(a[i].Uint64())
	}
	return mpcObj.TableLookup(aConv, tableId)
}

func (mpcObj *MPC) LessThanBits(a, b mpc_core.RMat) mpc_core.RVec {
	return mpcObj.LessThanBitsAux(a, b, false)
}

// b is public
func (mpcObj *MPC) LessThanBitsPublic(a, bpub mpc_core.RMat) mpc_core.RVec {
	return mpcObj.LessThanBitsAux(a, bpub, true)
}

// b is public if isPublic = true,
func (mpcObj *MPC) LessThanBitsAux(a, b mpc_core.RMat, isPublic bool) mpc_core.RVec {
	pid := mpcObj.Network.pid
	rtype := a.Type()
	n, k := a.Dims()
	if n != len(b) || k != len(b[0]) {
		panic("Input matrices a and b have mismatching dimensions")
	}

	// Compute XOR
	var x mpc_core.RMat
	if !isPublic {
		// Compute element-wise a * b
		x = mpcObj.SSMultElemMat(a, b)

		if pid > 0 {
			minusTwo := rtype.FromInt(-2)
			for i := range x {
				for j := range x[i] {
					x[i][j] = x[i][j].Mul(minusTwo).Add(a[i][j]).Add(b[i][j])
				}
			}
		}
	} else if pid > 0 {
		minusTwo := rtype.FromInt(-2)
		x = mpc_core.RMultElemMat(a, b)
		for i := range x {
			for j := range x[i] {
				x[i][j] = x[i][j].Mul(minusTwo).Add(a[i][j])
				if pid == 1 {
					x[i][j] = x[i][j].Add(b[i][j])
				}
			}
		}
	} else {
		x = mpc_core.InitRMat(rtype.Zero(), n, k)
	}

	f := mpcObj.PrefixOr(x)
	x = nil

	if pid > 0 {
		for i := 0; i < n; i++ {
			for j := k - 1; j >= 1; j-- {
				f[i][j] = f[i][j].Sub(f[i][j-1])
			}
		}
	}

	if isPublic {
		c := mpc_core.InitRVec(rtype.Zero(), n)
		if pid > 0 {
			for i := 0; i < n; i++ {
				for j := 0; j < k; j++ {
					c[i] = c[i].Add(f[i][j].Mul(b[i][j]))
				}
			}
		}
		return c
	}

	// Compute inner product for each pair of rows
	fr, fm := mpcObj.BeaverPartitionMat(f)
	br, bm := mpcObj.BeaverPartitionMat(b)
	c := mpcObj.BeaverMultElemMat(fr, fm, br, bm).Sum(0)
	c = mpcObj.BeaverReconstructVec(c)

	return c
}

func (mpcObj *MPC) NormalizerEvenExp2N(x mpc_core.RVec, k int) (mpc_core.RVec, mpc_core.RVec) {
	rtypeInput := x.Type().Zero()
	n := len(x)
	pid := mpcObj.GetPid()

	rtype2N := mpc_core.LElem2NBigIntZero
	rtypeBit := mpc_core.BElem(0)

	shareConversionBuffer := 15
	k += shareConversionBuffer

	if rtypeInput.TypeID() == mpc_core.LElem256UniqueID {
		rtype2N.SetModulusPowerOf2(256)
	} else if rtypeInput.TypeID() == mpc_core.LElem128UniqueID {
		rtype2N.SetModulusPowerOf2(127)
	} else {
		panic("Unsupported RElem type for NormalizerEvenExp2N")
	}

	// Pretend x is shared in the power-of-two ring
	// To cope with noise during conversion, left shift the value by shareConversionBuffer
	// and consider that to be the input. Increase k by the same amount to calculate the correct answer
	xConv2N := mpc_core.InitRVec(rtype2N, n)
	for i := range x {
		el := x[i]
		el = el.Mul(rtypeInput.FromInt(1 << shareConversionBuffer)) // shift = 2^(shareConversionBuffer)
		if rtypeInput.TypeID() == mpc_core.LElem256UniqueID {
			xConv2N[i] = xConv2N[i].FromBigInt(big.NewInt(0).Mod(el.(mpc_core.LElem256).ToBigInt(), el.Modulus()))
		} else if rtypeInput.TypeID() == mpc_core.LElem128UniqueID {
			xConv2N[i] = xConv2N[i].FromBigInt(big.NewInt(0).Mod(el.(mpc_core.LElem128).ToBigInt(), el.Modulus()))
		}
	}

	numBits := rtype2N.ModBitLength() - 1
	r, rBits := mpcObj.ShareRandomBits(rtype2N, rtypeBit, n, numBits, numBits, true)
	rBits = rBits.Transpose()

	// Reveal a = x + r
	r.Add(xConv2N)
	a := mpcObj.RevealSymVec(r)
	numBitsPacked := 1 + ((numBits - 1) / mpc_core.BElemBits)

	// Compute t
	tBits := mpc_core.InitRMat(rtypeBit.Zero(), numBitsPacked, n)
	tBits.Sub(rBits)
	if pid == mpcObj.GetHubPid() {
		tBits.AddScalar(rtypeBit.One())
	}

	// Get bit representation of a + 1
	a.AddScalar(rtype2N.FromInt(1))
	aBits := numToBits(rtypeBit, a, numBits)
	aBits = aBits.Transpose()

	// Bit-add aBits and tBits to obtain the bit representation of x
	xBits := mpcObj.BinaryAddPublic(aBits, tBits, numBits)

	// {
	// 	xBitsRev := mpcObj.RevealSymMat(xBits)
	// 	fmt.Printf(" xbits: %064b\n", xBitsRev[len(xBitsRev)-1][0].Uint64())
	// }

	// 0000101100 -> 1111000000
	pref := mpcObj.BinaryPrefixOr(xBits, numBits)

	// {
	// 	prefRev := mpcObj.RevealSymMat(pref)
	// 	fmt.Printf("  pref: %064b\n", prefRev[len(prefRev)-1][0].Uint64())
	// }

	numHalfBits := int(k / 2)

	r, rBits = mpcObj.ShareRandomBits(rtypeInput, rtypeBit, n*numHalfBits, 1, 1, true)

	numHalfElems := 1 + ((numHalfBits - 1) / mpc_core.BElemBits)
	halfBits := mpc_core.InitRMat(rtypeBit.Zero(), numHalfElems, n)

	// Extract half bits
	for i := 0; i < n; i++ {
		pos := 0
		for bitIndex := k - 2; bitIndex >= 0; bitIndex -= 2 {
			elemIndex := int(bitIndex / mpc_core.BElemBits)
			shift := uint64(bitIndex % mpc_core.BElemBits)
			bitMask := uint64(1) << shift
			bit := (pref[len(pref)-1-elemIndex][i].Uint64() & bitMask) >> shift
			halfBits[pos/mpc_core.BElemBits][i] = rtypeBit.FromUint64(halfBits[pos/mpc_core.BElemBits][i].Uint64() | (bit << (pos % mpc_core.BElemBits)))
			pos += 1
		}
	}

	// {
	// 	halfRev := mpcObj.RevealSymMat(halfBits)
	// 	fmt.Printf("  half: %064b\n", halfRev[len(halfRev)-1][0].Uint64())
	// }

	// Mask and reveal
	randInd := 0
	for i := 0; i < n; i++ {
		pos := 0
		for bitIndex := k - 2; bitIndex >= 0; bitIndex -= 2 {
			halfBits[pos/mpc_core.BElemBits][i] = rtypeBit.FromUint64(halfBits[pos/mpc_core.BElemBits][i].Uint64() ^ (rBits[randInd][0].Uint64() << (pos % mpc_core.BElemBits)))
			pos += 1
			randInd += 1
		}
	}
	halfBits = mpcObj.RevealSymMat(halfBits)

	// Weighted sum to get answers
	normalizer := mpc_core.InitRVec(rtypeInput.Zero(), n)
	normalizerSqrt := mpc_core.InitRVec(rtypeInput.Zero(), n)

	if pid == mpcObj.GetHubPid() {
		normalizer.AddScalar(rtypeInput.One())
		normalizerSqrt.AddScalar(rtypeInput.One())
	}

	coeff := mpc_core.InitRVec(rtypeInput, numHalfBits)
	coeffSqrt := mpc_core.InitRVec(rtypeInput, numHalfBits)
	prev, prevSqrt := rtypeInput.One(), rtypeInput.One()
	pos := 0
	for bitIndex := k - 2; bitIndex >= 0; bitIndex -= 2 { // Precompute coefficients
		cur := prev.Mul(rtypeInput.FromInt(4))
		curSqrt := prevSqrt.Mul(rtypeInput.FromInt(2))
		coeff[pos] = cur.Sub(prev)
		coeffSqrt[pos] = curSqrt.Sub(prevSqrt)
		prev, prevSqrt = cur, curSqrt
		pos += 1
	}

	randInd = 0
	for i := 0; i < n; i++ {
		pos := 0
		for bitIndex := k - 2; bitIndex >= 0; bitIndex -= 2 {
			flip := halfBits[pos/mpc_core.BElemBits][i].Uint64() & (uint64(1) << (pos % mpc_core.BElemBits))
			randShare := r[randInd]
			if flip > 0 { // Take 1 - [r], which gives the secret bit
				if pid == mpcObj.GetHubPid() {
					randShare = randShare.One().Sub(randShare)
				} else {
					randShare = randShare.Neg()
				}
			}

			normalizer[i] = normalizer[i].Add(randShare.Mul(coeff[pos]))
			normalizerSqrt[i] = normalizerSqrt[i].Add(randShare.Mul(coeffSqrt[pos]))

			pos += 1
			randInd += 1
		}
	}

	return normalizer, normalizerSqrt
}

// Returns power-of-two multipliers b and sqrt(b) (for each a)
// where the bit length of a * b becomes k or k-1.
// Requires input a to be strictly positive and k to be even
// For statistical security, also requires ring modulus to be at least
// 30 bits longer than the maximum data bit length k
func (mpcObj *MPC) NormalizerEvenExp(a mpc_core.RVec, k int) (mpc_core.RVec, mpc_core.RVec) {
	pid := mpcObj.GetPid()
	n := len(a)
	rtype := a.Type()
	rtypeBit := mpc_core.SElemDS(0)

	r, rBits := mpcObj.ShareRandomBits(rtype, rtypeBit, n, k, k+30, true)

	e := a.Copy()
	e.Add(r)
	r = nil

	e = mpcObj.RevealSymVec(e) // Reveal a + r
	eBits := numToBits(rtypeBit, e, k)
	e = nil

	c := mpcObj.LessThanBitsPublic(rBits, eBits)
	if pid > 0 { // Compute 1 - c
		var cNew mpc_core.RVec
		if pid == 1 {
			cNew = mpc_core.InitRVec(rtypeBit.One(), n)
		} else {
			cNew = mpc_core.InitRVec(rtypeBit.Zero(), n)
		}
		cNew.Sub(c)
		c = cNew
	}

	ep := mpc_core.InitRMat(rtypeBit.Zero(), n, k+1)
	two := rtypeBit.FromInt(2)
	if pid > 0 {
		for i := 0; i < n; i++ {
			ep[i][0] = c[i]
			for j := 1; j < k+1; j++ {
				ep[i][j] = rtypeBit.One().Sub(eBits[i][j-1].Mul(two)).Mul(rBits[i][j-1])
				if pid == 1 {
					ep[i][j] = ep[i][j].Add(eBits[i][j-1])
				}
			}
		}
	}
	c = nil

	E := mpcObj.PrefixOr(ep)

	tpNeg := mpc_core.InitRMat(rtypeBit.Zero(), n, k)
	if pid > 0 {
		for i := 0; i < n; i++ {
			for j := 0; j < k; j++ {
				tpNeg[i][j] = E[i][j].Sub(rBits[i][j].Mul(rtypeBit.One().Sub(eBits[i][j])))
			}
		}
	}
	E = nil

	TNeg := mpcObj.PrefixOr(tpNeg)
	tpNeg = nil

	halfLen := int(k / 2)
	efir := mpc_core.RMultElemMat(eBits, TNeg)
	rfir := mpcObj.SSMultElemMat(rBits, TNeg)
	eBits, rBits = nil, nil

	doubleFlag := mpcObj.LessThanBits(efir, rfir)
	efir, rfir = nil, nil

	oddBits := mpc_core.InitRMat(rtypeBit.Zero(), n, halfLen)
	evenBits := mpc_core.InitRMat(rtypeBit.Zero(), n, halfLen)
	if pid > 0 {
		for i := 0; i < n; i++ {
			for j := 0; j < halfLen; j++ {
				oddBits[i][j] = oddBits[i][j].Sub(TNeg[i][2*j+1])
				if pid == 1 {
					oddBits[i][j] = oddBits[i][j].Add(rtypeBit.One())
				}
				if 2*j+2 < k {
					evenBits[i][j] = evenBits[i][j].Sub(TNeg[i][2*j+2])
					if pid == 1 {
						evenBits[i][j] = evenBits[i][j].Add(rtypeBit.One())
					}
				} // zero otherwise
			}
		}
	}
	TNeg = nil

	oddBitSum := oddBits.Sum(0)
	evenBitSum := evenBits.Sum(0)
	if pid == 1 {
		oddBitSum.AddScalar(rtypeBit.One())
		evenBitSum.AddScalar(rtypeBit.One())
	}
	oddBits, evenBits = nil, nil

	// If doubleFlag = true, then use oddBits, otherwise use evenBits
	diff := oddBitSum.Copy()
	diff.Sub(evenBitSum)
	diff = mpcObj.SSMultElemVec(doubleFlag, diff)
	doubleFlag = nil

	chosenBitSum := evenBitSum.Copy()
	chosenBitSum.Add(diff)
	oddBitSum, evenBitSum, diff = nil, nil, nil

	bMat := mpcObj.TableLookupWithShareConversion(rtype, chosenBitSum, 1)

	return bMat[1], bMat[0]
}

func (mpcObj *MPC) BinaryPrefixOr(a mpc_core.RMat, numBits int) mpc_core.RMat {
	rtype := a.Type()
	pid := mpcObj.GetPid()

	if rtype.TypeID() != mpc_core.BElemUniqueID {
		panic("BinaryPrefixOr defined only for BElem")
	}

	// Mask leading bits (set secret bits to 0)
	numBitsHead := numBits % mpc_core.BElemBits
	if numBitsHead > 0 {
		tipMask := ^uint64(0) << numBitsHead // 11..1100..00
		for i := range a[0] {
			a[0][i] = mpc_core.BElem(a[0][i].Uint64() & ^tipMask)
		}
	}

	// Invert all secret bits: 00001010011 -> 11110101100
	out := a.Copy()
	if pid == mpcObj.GetHubPid() {
		for i := range a {
			for j := range a[i] {
				out[i][j] = mpc_core.BElem(^a[i][j].Uint64())
			}
		}
	}

	if mpcObj.GetPid() > 0 {
		// chosen := 0
		// fmt.Printf("dpro: %064b (round=%d) input\n", mpcObj.RevealSym(out[0][chosen]).Uint64(), 0)
		// fmt.Printf("dpro: %064b%064b%064b%064b (round=%d) input\n", mpcObj.RevealSym(out[0][chosen]).Uint64(), mpcObj.RevealSym(out[1][chosen]).Uint64(), mpcObj.RevealSym(out[2][chosen]).Uint64(), mpcObj.RevealSym(out[3][chosen]).Uint64(), 0)
	}

	numRounds := int(math.Ceil(math.Log2(float64(numBits))))
	for round := 1; round <= numRounds; round++ {
		if round <= 6 { // Shift stays within each uint64

			mask := (^uint64(0) / ((1 << (1 << (round))) - 1)) << ((1 << (round - 1)) - 1)
			maskMove := mask << 1
			for sh := 0; sh < round-1; sh++ {
				maskMove |= (maskMove << (1 << sh))
			}
			mask <<= 1
			maskMove = ^maskMove

			move := out.Copy()

			// Pick product bits and move to right with replication
			// E.g.: 0001000000010000 -> 0000111100001111
			for i := range move {
				for j := range move[i] {
					val := (move[i][j].Uint64() & mask) >> 1
					for sh := 0; sh < round-1; sh++ {
						val |= (val >> (1 << sh))
					}
					move[i][j] = rtype.FromUint64(val)
				}
			}

			// Multiply with previous carry bits
			sr, sm := mpcObj.BeaverPartitionMat(move)
			or, om := mpcObj.BeaverPartitionMat(out)

			mult := mpcObj.BeaverMultElemMat(or, om, sr, sm)

			mult = mpcObj.BeaverReconstructMat(mult)

			// Update carry bits
			for i := range out {
				for j := range out[i] {
					out[i][j] = rtype.FromUint64((out[i][j].Uint64() & ^maskMove) | (mult[i][j].Uint64() & maskMove))
				}
			}

		} else { // Need to move across uint64

			numElems := len(out)
			n := len(out[0])

			mod := 1 << (round - 6) // 2, 4, 8, 16, ...
			chosen := (mod / 2) - 1 // 0, 1, 3, 7, ...

			numIn, numOut := 0, 0
			for j := 0; j < numElems; j++ {
				if j%mod == chosen {
					numIn++
				} else if j%mod > chosen {
					numOut++
				}
			}

			propIn := make(mpc_core.RMat, numIn)
			propOut := make(mpc_core.RMat, numOut)

			indIn, indOut := 0, 0
			for j := 0; j < numElems; j++ {
				if j%mod == chosen {
					propIn[indIn] = out[j].Copy()
					for s := 0; s < n; s++ {
						maskHead := uint64(1)
						tmp := propIn[indIn][s].Uint64() & maskHead
						for sh := 0; sh < 6; sh++ {
							tmp |= (tmp << (1 << sh))
						}
						propIn[indIn][s] = rtype.FromUint64(tmp)
					}
					indIn++
				} else if j%mod > chosen {
					propOut[indOut] = out[j].Copy()
					indOut++
				}
			}

			dPSr, dPSm := mpcObj.BeaverPartitionMat(propIn)
			dPr, dPm := mpcObj.BeaverPartitionMat(propOut)

			propOut.Clear()

			indIn, indOut = -1, 0
			for j := 0; j < numElems; j++ {
				if j%mod == chosen {
					indIn++
				} else if j%mod > chosen {
					propOut[indOut] = mpcObj.BeaverMultElemVec(dPSr[indIn], dPSm[indIn], dPr[indOut], dPm[indOut])
					indOut++
				}
			}

			propOut = mpcObj.BeaverReconstructMat(propOut)

			indOut = 0
			for j := 0; j < numElems; j++ {
				if j%mod > chosen {
					out[j] = propOut[indOut].Copy()
					indOut++
				}
			}
		}

		if mpcObj.GetPid() > 0 {
			// chosen := 0
			// fmt.Printf("dpro: %064b (round=%d) input\n", mpcObj.RevealSym(out[0][chosen]).Uint64(), round)
			// fmt.Printf("dpro: %064b%064b%064b%064b (round=%d) input\n", mpcObj.RevealSym(out[0][chosen]).Uint64(), mpcObj.RevealSym(out[1][chosen]).Uint64(), mpcObj.RevealSym(out[2][chosen]).Uint64(), mpcObj.RevealSym(out[3][chosen]).Uint64(), round)
		}
	}

	// if mpcObj.GetPid() > 0 {
	// 	chosen := 0
	// 	fmt.Printf("  dpro: %064b (final)\n", mpcObj.RevealSym(dProp[0][chosen]).Uint64())
	// 	fmt.Printf("  dgen: %064b (final)\n", mpcObj.RevealSym(dGen[0][chosen]).Uint64())
	// 	// fmt.Printf("dpro: %064b%064b%064b%064b (numBits=%d) input\n", mpcObj.RevealSym(dProp[0][chosen]).Uint64(), mpcObj.RevealSym(dProp[1][chosen]).Uint64(), mpcObj.RevealSym(dProp[2][chosen]).Uint64(), mpcObj.RevealSym(dProp[3][chosen]).Uint64(), numBits)
	// 	// fmt.Printf("dgen: %064b%064b%064b%064b (numBits=%d) input\n", mpcObj.RevealSym(dGen[0][chosen]).Uint64(), mpcObj.RevealSym(dGen[1][chosen]).Uint64(), mpcObj.RevealSym(dGen[2][chosen]).Uint64(), mpcObj.RevealSym(dGen[3][chosen]).Uint64(), numBits)
	// }

	return out
}

func (mpcObj *MPC) BinaryAddPublic(aPub mpc_core.RMat, b mpc_core.RMat, numBits int) mpc_core.RMat {
	pid := mpcObj.GetPid()
	carry := mpcObj.CarryOverPublic(aPub, b, numBits)
	res := b.Copy()
	res.Add(carry)
	if pid == mpcObj.GetHubPid() {
		res.Add(aPub)
	}
	return res
}

// Each column of aPub and b encodes a bitwise-shared number (mpc_core.BElem) with length numBits
// Outputs all carry bits of the binary addition of each pair of columns (in each output column)
// Reference: https://faui1-files.cs.fau.de/filepool/publications/octavian_securescm/SecureSCM-D.9.2.pdf
func (mpcObj *MPC) CarryOverPublic(aPub mpc_core.RMat, b mpc_core.RMat, numBits int) mpc_core.RMat {
	// Carry propagation bits
	dProp := b.Copy()
	if mpcObj.GetPid() == mpcObj.GetHubPid() {
		dProp.Add(aPub)
	}

	// Carry generation bits
	dGen := b.Copy()
	dGen.MulElem(aPub)

	return mpcObj.SuffixCarryAux(dProp, dGen, numBits)
}

func (mpcObj *MPC) SuffixCarryAux(dProp mpc_core.RMat, dGen mpc_core.RMat, numBits int) mpc_core.RMat {
	rtype := dProp.Type()

	// Mask leading bits (set secret bits to 1 for dProp and 0 for dGen)
	numBitsHead := numBits % mpc_core.BElemBits
	if numBitsHead > 0 {
		tipMask := ^uint64(0) << numBitsHead // 11..1100..00
		for i := range dProp[0] {
			if mpcObj.GetPid() == mpcObj.GetHubPid() {
				dProp[0][i] = mpc_core.BElem(dProp[0][i].Uint64() | tipMask)
			} else {
				dProp[0][i] = mpc_core.BElem(dProp[0][i].Uint64() & ^tipMask)
			}
			dGen[0][i] = mpc_core.BElem(dGen[0][i].Uint64() & ^tipMask)
		}
	}

	// if mpcObj.GetPid() > 0 {
	// 	chosen := 0
	// 	fmt.Printf("  dpro: %064b (round=%d)\n", mpcObj.RevealSym(dProp[0][chosen]).Uint64(), 0)
	// 	fmt.Printf("  dgen: %064b (round=%d)\n", mpcObj.RevealSym(dGen[0][chosen]).Uint64(), 0)
	// 	// fmt.Printf("dpro: %064b%064b%064b%064b (numBits=%d) input\n", mpcObj.RevealSym(dProp[0][chosen]).Uint64(), mpcObj.RevealSym(dProp[1][chosen]).Uint64(), mpcObj.RevealSym(dProp[2][chosen]).Uint64(), mpcObj.RevealSym(dProp[3][chosen]).Uint64(), numBits)
	// 	// fmt.Printf("dgen: %064b%064b%064b%064b (numBits=%d) input\n", mpcObj.RevealSym(dGen[0][chosen]).Uint64(), mpcObj.RevealSym(dGen[1][chosen]).Uint64(), mpcObj.RevealSym(dGen[2][chosen]).Uint64(), mpcObj.RevealSym(dGen[3][chosen]).Uint64(), numBits)
	// }

	numRounds := int(math.Ceil(math.Log2(float64(numBits))))
	for round := 1; round <= numRounds; round++ {
		if round <= 6 { // Shift stays within each uint64

			mask := (^uint64(0) / ((1 << (1 << (round))) - 1)) << ((1 << (round - 1)) - 1)
			maskMove := mask << 1
			for sh := 0; sh < round-1; sh++ {
				maskMove |= (maskMove << (1 << sh))
			}

			dPropMove := dProp.Copy()
			dGenMove := dGen.Copy()

			// Pick product bits and move to left with replication
			// E.g.: 0000100000001000 -> 1111000011110000
			for i := range dPropMove {
				for j := range dPropMove[i] {
					propMove := (dPropMove[i][j].Uint64() & mask) << 1
					genMove := (dGenMove[i][j].Uint64() & mask) << 1
					for sh := 0; sh < round-1; sh++ {
						propMove |= (propMove << (1 << sh))
						genMove |= (genMove << (1 << sh))
					}
					dPropMove[i][j] = rtype.FromUint64(propMove)
					dGenMove[i][j] = rtype.FromUint64(genMove)
				}
			}

			// Multiply with previous carry bits
			dPSr, dPSm := mpcObj.BeaverPartitionMat(dPropMove)
			dGSr, dGSm := mpcObj.BeaverPartitionMat(dGenMove)
			dPr, dPm := mpcObj.BeaverPartitionMat(dProp)

			dPropOut := mpcObj.BeaverMultElemMat(dPSr, dPSm, dPr, dPm)
			dGenOut := mpcObj.BeaverMultElemMat(dGSr, dGSm, dPr, dPm)

			dPropOut = mpcObj.BeaverReconstructMat(dPropOut)
			dGenOut = mpcObj.BeaverReconstructMat(dGenOut)
			dGenOut.Add(dGen)

			// Update carry bits
			for i := range dPropMove {
				for j := range dPropMove[i] {
					dProp[i][j] = rtype.FromUint64((dProp[i][j].Uint64() & ^maskMove) | (dPropOut[i][j].Uint64() & maskMove))
					dGen[i][j] = rtype.FromUint64((dGen[i][j].Uint64() & ^maskMove) | (dGenOut[i][j].Uint64() & maskMove))
				}
			}

		} else { // Need to move across uint64

			numElems := len(dProp)
			n := len(dProp[0])

			mod := 1 << (round - 6) // 2, 4, 8, 16, ...
			chosen := (mod / 2) - 1 // 0, 1, 3, 7, ...

			numIn, numOut := 0, 0
			for j := 0; j < numElems; j++ {
				if j%mod == chosen {
					numIn++
				} else if j%mod > chosen {
					numOut++
				}
			}

			propIn := make(mpc_core.RMat, numIn)
			genIn := make(mpc_core.RMat, numIn)
			propOut := make(mpc_core.RMat, numOut)
			genOut := make(mpc_core.RMat, numOut)

			indIn, indOut := 0, 0
			for j := 0; j < numElems; j++ {
				if j%mod == chosen {
					propIn[indIn] = dProp[numElems-1-j].Copy()
					genIn[indIn] = dGen[numElems-1-j].Copy()
					for s := 0; s < n; s++ {
						maskHead := uint64(1) << 63
						tmp1 := propIn[indIn][s].Uint64() & maskHead
						tmp2 := genIn[indIn][s].Uint64() & maskHead
						for sh := 0; sh < 6; sh++ {
							tmp1 |= (tmp1 >> (1 << sh))
							tmp2 |= (tmp2 >> (1 << sh))
						}
						propIn[indIn][s] = rtype.FromUint64(tmp1)
						genIn[indIn][s] = rtype.FromUint64(tmp2)
					}
					indIn++
				} else if j%mod > chosen {
					propOut[indOut] = dProp[numElems-1-j].Copy()
					genOut[indOut] = dGen[numElems-1-j].Copy()
					indOut++
				}
			}

			dPSr, dPSm := mpcObj.BeaverPartitionMat(propIn)
			dGSr, dGSm := mpcObj.BeaverPartitionMat(genIn)
			dPr, dPm := mpcObj.BeaverPartitionMat(propOut)

			propOut.Clear()
			genOut.Clear()

			indIn, indOut = -1, 0
			for j := 0; j < numElems; j++ {
				if j%mod == chosen {
					indIn++
				} else if j%mod > chosen {
					propOut[indOut] = mpcObj.BeaverMultElemVec(dPSr[indIn], dPSm[indIn], dPr[indOut], dPm[indOut])
					genOut[indOut] = mpcObj.BeaverMultElemVec(dGSr[indIn], dGSm[indIn], dPr[indOut], dPm[indOut])
					indOut++
				}
			}

			propOut = mpcObj.BeaverReconstructMat(propOut)
			genOut = mpcObj.BeaverReconstructMat(genOut)

			indOut = 0
			for j := 0; j < numElems; j++ {
				if j%mod > chosen {
					dProp[numElems-1-j] = propOut[indOut].Copy()
					dGen[numElems-1-j].Add(genOut[indOut])
					indOut++
				}
			}
		}

		// if mpcObj.GetPid() > 0 {
		// 	chosen := 0
		// 	fmt.Printf("  dpro: %064b (round=%d)\n", mpcObj.RevealSym(dProp[0][chosen]).Uint64(), round)
		// 	fmt.Printf("  dgen: %064b (round=%d)\n", mpcObj.RevealSym(dGen[0][chosen]).Uint64(), round)
		// 	// fmt.Printf("dpro: %064b%064b%064b%064b (numBits=%d) input\n", mpcObj.RevealSym(dProp[0][chosen]).Uint64(), mpcObj.RevealSym(dProp[1][chosen]).Uint64(), mpcObj.RevealSym(dProp[2][chosen]).Uint64(), mpcObj.RevealSym(dProp[3][chosen]).Uint64(), numBits)
		// 	// fmt.Printf("dgen: %064b%064b%064b%064b (numBits=%d) input\n", mpcObj.RevealSym(dGen[0][chosen]).Uint64(), mpcObj.RevealSym(dGen[1][chosen]).Uint64(), mpcObj.RevealSym(dGen[2][chosen]).Uint64(), mpcObj.RevealSym(dGen[3][chosen]).Uint64(), numBits)
		// }
	}

	if numBitsHead > 0 {
		tipMask := ^uint64(0) << numBitsHead // 11..1100..00
		for i := range dGen[0] {
			dGen[0][i] = rtype.FromUint64((dGen[0][i].Uint64() & ^tipMask))
		}
	}

	// Left shift by one to align the carry bits with target positions
	firstBit := make([]bool, len(dGen[0]))
	headMask := uint64(1) << 63
	for j := len(dGen) - 1; j >= 0; j-- {
		for i := range dGen[j] {
			val := dGen[j][i].Uint64()
			head := val&headMask > 0
			val <<= 1
			if j < len(dGen)-1 && firstBit[i] { // potential carry from previous uint64
				val += 1
			}
			dGen[j][i] = rtype.FromUint64(val)
			firstBit[i] = head
		}
	}

	// if mpcObj.GetPid() > 0 {
	// 	chosen := 0
	// 	fmt.Printf("  dpro: %064b (final)\n", mpcObj.RevealSym(dProp[0][chosen]).Uint64())
	// 	fmt.Printf("  dgen: %064b (final)\n", mpcObj.RevealSym(dGen[0][chosen]).Uint64())
	// 	// fmt.Printf("dpro: %064b%064b%064b%064b (numBits=%d) input\n", mpcObj.RevealSym(dProp[0][chosen]).Uint64(), mpcObj.RevealSym(dProp[1][chosen]).Uint64(), mpcObj.RevealSym(dProp[2][chosen]).Uint64(), mpcObj.RevealSym(dProp[3][chosen]).Uint64(), numBits)
	// 	// fmt.Printf("dgen: %064b%064b%064b%064b (numBits=%d) input\n", mpcObj.RevealSym(dGen[0][chosen]).Uint64(), mpcObj.RevealSym(dGen[1][chosen]).Uint64(), mpcObj.RevealSym(dGen[2][chosen]).Uint64(), mpcObj.RevealSym(dGen[3][chosen]).Uint64(), numBits)
	// }

	return dGen
}

// Each column of aPub and b encodes a bitwise-shared number (mpc_core.BElem) with length numBits
// Outputs the final carry bit of the binary addition of each pair of columns (as LSB of mpc_core.BElem)
// Reference: https://faui1-files.cs.fau.de/filepool/publications/octavian_securescm/SecureSCM-D.9.2.pdf
func (mpcObj *MPC) CarryOutPublic(aPub mpc_core.RMat, b mpc_core.RMat, numBits int) mpc_core.RVec {
	// Carry propagation bits
	dProp := b.Copy()
	if mpcObj.GetPid() == mpcObj.GetHubPid() {
		dProp.Add(aPub)
	}

	// Carry generation bits
	dGen := b.Copy()
	dGen.MulElem(aPub)

	return mpcObj.CarryOutAux(dProp, dGen, numBits)
}

func (mpcObj *MPC) CarryOutAux(dProp mpc_core.RMat, dGen mpc_core.RMat, numBits int) mpc_core.RVec {
	// log.LLvl1("CarryOutAux", "pid", mpcObj.GetPid(), "dProp", len(dProp), len(dProp[0]), "dGen", len(dGen), len(dGen[0]))

	// chosen := 2
	// if mpcObj.GetPid() > 0 {
	// 	fmt.Printf("dpro: %064b%064b%064b%064b (numBits=%d) input\n", mpcObj.RevealSym(dProp[0][chosen]).Uint64(), mpcObj.RevealSym(dProp[1][chosen]).Uint64(), mpcObj.RevealSym(dProp[2][chosen]).Uint64(), mpcObj.RevealSym(dProp[3][chosen]).Uint64(), numBits)
	// 	fmt.Printf("dgen: %064b%064b%064b%064b (numBits=%d) input\n", mpcObj.RevealSym(dGen[0][chosen]).Uint64(), mpcObj.RevealSym(dGen[1][chosen]).Uint64(), mpcObj.RevealSym(dGen[2][chosen]).Uint64(), mpcObj.RevealSym(dGen[3][chosen]).Uint64(), numBits)
	// }

	if numBits == 1 {
		rv := dGen[0]
		for i := range rv {
			rv[i] = mpc_core.BElem(rv[i].Uint64() & 1)
		}
		return rv
	} else {
		// Mask leading bits (set secret bits to 1 for dProp and 0 for dGen)
		numBitsHead := numBits % mpc_core.BElemBits
		if numBitsHead > 0 {
			tipMask := ^uint64(0) << numBitsHead // 11..1100..00
			for i := range dProp[0] {
				if mpcObj.GetPid() == mpcObj.GetHubPid() {
					dProp[0][i] = mpc_core.BElem(dProp[0][i].Uint64() | tipMask)
				} else {
					dProp[0][i] = mpc_core.BElem(dProp[0][i].Uint64() & ^tipMask)
				}
				dGen[0][i] = mpc_core.BElem(dGen[0][i].Uint64() & ^tipMask)
			}
		}

		// if mpcObj.GetPid() > 0 {
		// 	fmt.Printf("dpro: %064b (numBits=%d) masked\n", mpcObj.RevealSym(dProp[0][chosen]).Uint64(), numBits)
		// 	fmt.Printf("dgen: %064b (numBits=%d) masked\n", mpcObj.RevealSym(dGen[0][chosen]).Uint64(), numBits)
		// }

		dPropShift := dProp.Copy()
		dGenShift := dGen.Copy()
		for i := range dProp {
			for j := range dProp[i] {
				dPropShift[i][j] = mpc_core.BElem(dPropShift[i][j].Uint64() >> 1)
				dGenShift[i][j] = mpc_core.BElem(dGenShift[i][j].Uint64() >> 1)
			}
		}

		dPSr, dPSm := mpcObj.BeaverPartitionMat(dPropShift)
		dPr, dPm := mpcObj.BeaverPartitionMat(dProp)
		dGr, dGm := mpcObj.BeaverPartitionMat(dGen)

		dPropOut := mpcObj.BeaverMultElemMat(dPSr, dPSm, dPr, dPm)
		dGenOut := mpcObj.BeaverMultElemMat(dPSr, dPSm, dGr, dGm)

		dPropOut = mpcObj.BeaverReconstructMat(dPropOut)
		dGenOut = mpcObj.BeaverReconstructMat(dGenOut)
		dGenOut.Add(dGenShift)

		// if mpcObj.GetPid() > 0 {
		// 	fmt.Printf("------%064b\n", (^uint64(0))/3)
		// 	fmt.Printf("dpro: %064b (numBits=%d) pair merge\n", mpcObj.RevealSym(dPropOut[0][chosen]).Uint64(), numBits)
		// 	fmt.Printf("dgen: %064b (numBits=%d) pair merge\n", mpcObj.RevealSym(dGenOut[0][chosen]).Uint64(), numBits)
		// }

		// 01010101 -> 00001111 (Collect all odd position bits in the second half of the integer)
		for i := 1; i < int(math.Log2(float64(mpc_core.BElemBits))); i++ {
			mask1 := (^uint64(0) / ((uint64(1) << (1 << (i + 1))) - 1)) * (((uint64(1) << (1 << i)) - 1) << (1 << i))
			mask2 := ^((^uint64(0) / ((uint64(1) << (1 << i)) - 1)) * (((uint64(1) << (1 << (i - 1))) - 1) << (1 << (i - 1))))
			maskMove := mask1 & mask2
			maskKeep := ^mask1 & mask2
			// if mpcObj.GetPid() > 0 {
			// 	fmt.Printf("-msk1-%064b (i=%d)\n", maskMove, i)
			// 	fmt.Printf("-msk2-%064b (i=%d)\n", maskKeep, i)
			// }

			for j := range dPropOut {
				for k := range dPropOut[j] {
					dPropOut[j][k] = mpc_core.BElem(((maskMove & dPropOut[j][k].Uint64()) >> (1 << (i - 1))) | (maskKeep & dPropOut[j][k].Uint64()))
					dGenOut[j][k] = mpc_core.BElem(((maskMove & dGenOut[j][k].Uint64()) >> (1 << (i - 1))) | (maskKeep & dGenOut[j][k].Uint64()))
				}
			}
		}

		// if mpcObj.GetPid() > 0 {
		// 	fmt.Printf("dpro: %064b (numBits=%d) collect\n", mpcObj.RevealSym(dPropOut[0][chosen]).Uint64(), numBits)
		// 	fmt.Printf("dgen: %064b (numBits=%d) collect\n", mpcObj.RevealSym(dGenOut[0][chosen]).Uint64(), numBits)
		// }

		// Pack adjacent pairs of integers into a single integer: (00001111, 00001111) -> (00000000, 11111111)
		halfBits := uint64(mpc_core.BElemBits) / 2
		halfMask := (uint64(1) << halfBits) - 1

		numBitsOut := (numBits + 1) / 2
		nelemOut := 1 + ((numBitsOut - 1) / mpc_core.BElemBits)
		nelemIn := 1 + ((numBits - 1) / mpc_core.BElemBits)

		startIndex := 0
		outIndex := 0
		dPropOutHalved := make(mpc_core.RMat, nelemOut)
		dGenOutHalved := make(mpc_core.RMat, nelemOut)
		if nelemIn%2 != 0 {
			startIndex += 1
			outIndex += 1
			dPropOutHalved[0] = dPropOut[0]
			dGenOutHalved[0] = dGenOut[0]
		}
		for i := startIndex; i < nelemIn; i += 2 {
			dPropOutHalved[outIndex] = dPropOut[i+1].Copy()
			dGenOutHalved[outIndex] = dGenOut[i+1].Copy()
			for k := range dPropOut[i] {
				dPropOutHalved[outIndex][k] = mpc_core.BElem(((dPropOut[i][k].Uint64() & halfMask) << halfBits) | (dPropOut[i+1][k].Uint64() & halfMask))
				dGenOutHalved[outIndex][k] = mpc_core.BElem(((dGenOut[i][k].Uint64() & halfMask) << halfBits) | (dGenOut[i+1][k].Uint64() & halfMask))
			}
			outIndex += 1
		}

		return mpcObj.CarryOutAux(dPropOutHalved, dGenOutHalved, numBitsOut)
	}
}

// Assumes x is LElem256 or LElem128 for now
func (mpcObj *MPC) IsPositive2N(x mpc_core.RVec) mpc_core.RVec {
	rtypeInput := x.Type().Zero()
	n := len(x)
	pid := mpcObj.GetPid()

	rtype2N := mpc_core.LElem2NBigIntZero
	rtypeBit := mpc_core.BElem(0)

	if rtypeInput.TypeID() == mpc_core.LElem256UniqueID {
		rtype2N.SetModulusPowerOf2(256)
	} else if rtypeInput.TypeID() == mpc_core.LElem128UniqueID {
		rtype2N.SetModulusPowerOf2(127)
	} else {
		panic("Unsupported RElem type for IsPositive2N")
	}

	// Pretend x is shared in the power-of-two ring
	// Transformation: pretend(shift * (x + 1)) - 1, where shift > (number of parties) * (modulus gap)
	// This preserves the sign of x even with the conversion noise
	xConv2N := mpc_core.InitRVec(rtype2N, n)
	for i := range x {
		el := x[i]
		if pid == mpcObj.GetHubPid() {
			el = el.Add(rtypeInput.One())
		}
		el = el.Mul(rtypeInput.FromInt(1 << 15)) // shift = 2^15
		if rtypeInput.TypeID() == mpc_core.LElem256UniqueID {
			xConv2N[i] = xConv2N[i].FromBigInt(big.NewInt(0).Mod(el.(mpc_core.LElem256).ToBigInt(), el.Modulus()))
		} else if rtypeInput.TypeID() == mpc_core.LElem128UniqueID {
			xConv2N[i] = xConv2N[i].FromBigInt(big.NewInt(0).Mod(el.(mpc_core.LElem128).ToBigInt(), el.Modulus()))
		}
		if pid == mpcObj.GetHubPid() {
			xConv2N[i] = xConv2N[i].Sub(rtype2N.One())
		}
	}

	// {
	// 	if pid > 0 {
	// 		fmt.Println("## X CONV: ", mpcObj.RevealSymVec(x).ToFloat(20)[:5])
	// 		fmt.Println("## X CONV: ", mpcObj.RevealSymVec(xConv2N).ToFloat(20)[:5])
	// 	}
	// }
	// panic("STOP")

	numBits := rtype2N.ModBitLength() - 1
	r, rBits := mpcObj.ShareRandomBits(rtype2N, rtypeBit, n, numBits, numBits, true)
	rBits = rBits.Transpose()

	// {
	// 	if pid > 0 {
	// 		fmt.Printf("## r     : %0256b\n", mpcObj.RevealSymVec(r)[0].(mpc_core.LElem2NBigInt).ToBigInt())
	// 		fmt.Printf("## r bits: %064b%064b%064b%064b\n", mpcObj.RevealSymMat(rBits).Transpose()[0][0].Uint64(), mpcObj.RevealSymMat(rBits).Transpose()[0][1].Uint64(), mpcObj.RevealSymMat(rBits).Transpose()[0][2].Uint64(), mpcObj.RevealSymMat(rBits).Transpose()[0][3].Uint64())
	// 	}
	// }
	// panic("STOP")

	// chosen := 2
	// x0 := mpcObj.RevealSym(xConv2N[chosen]).Uint64()
	// r0 := mpcObj.RevealSym(r[0]).Uint64()
	// {
	// 	if pid > 0 {
	// 		fmt.Printf("x[0]: %064b\n", x0)
	// 	}
	// }

	// Reveal a = x + r
	r.Add(xConv2N)
	a := mpcObj.RevealSymVec(r)
	numBitsPacked := 1 + ((numBits - 1) / mpc_core.BElemBits)

	// {
	// 	if pid > 0 {
	// 		fmt.Printf("a[0]: %064b\n", a[chosen].Uint64())
	// 	}
	// }

	// Compute t
	tBits := mpc_core.InitRMat(rtypeBit.Zero(), numBitsPacked, n)
	tBits.Sub(rBits)
	if pid == mpcObj.GetHubPid() {
		tBits.AddScalar(rtypeBit.One())
	}

	// Get bit representation of a + 1
	a.AddScalar(rtype2N.FromInt(1))
	aBits := numToBits(rtypeBit, a, numBits)
	aBits = aBits.Transpose()

	// a0 := aBits[0][chosen].Uint64()
	// t0 := mpcObj.RevealSym(tBits[0][chosen]).Uint64()
	// {
	// 	if pid > 0 {
	// 		fmt.Printf("rbit: %064b%064b%064b%064b\n", mpcObj.RevealSym(rBits[0][chosen]).Uint64(), mpcObj.RevealSym(rBits[1][chosen]).Uint64(), mpcObj.RevealSym(rBits[2][chosen]).Uint64(), mpcObj.RevealSym(rBits[3][chosen]).Uint64())
	// 		fmt.Printf("tbit: %064b%064b%064b%064b\n", mpcObj.RevealSym(tBits[0][chosen]).Uint64(), mpcObj.RevealSym(tBits[1][chosen]).Uint64(), mpcObj.RevealSym(tBits[2][chosen]).Uint64(), mpcObj.RevealSym(tBits[3][chosen]).Uint64())
	// 		fmt.Printf("abit: %064b%064b%064b%064b\n", aBits[0][chosen].Uint64(), aBits[1][chosen].Uint64(), aBits[2][chosen].Uint64(), aBits[3][chosen].Uint64())
	// 	}
	// }

	// Bit-add aBits and tBits to obtain the incoming carry bit at the MSB position
	carryMSB := mpcObj.CarryOutPublic(aBits, tBits, numBits-1)
	// {
	// 	if pid > 0 {
	// 		fmt.Printf("===carryMSB: %d\n a+t: %064b", mpcObj.RevealSym(carryMSB[chosen]).Uint64(), t0+a0)
	// 		fmt.Printf("Result: %d, Xor: %d\n", ((t0+a0)&(1<<numBits-1)>>(numBits-1))%2, (((t0^a0)&(1<<numBits-1))>>(numBits-1))%2)
	// 	}
	// }

	msbPos := (numBits - 1) % mpc_core.BElemBits
	aMSB := mpc_core.InitRVec(rtypeBit.Zero(), n)
	tMSB := mpc_core.InitRVec(rtypeBit.Zero(), n)
	MSB := mpc_core.InitRVec(rtypeBit.Zero(), n)
	for i := 0; i < n; i++ {
		aMSB[i] = rtypeBit.FromUint64(uint64(aBits[0][i].GetBit(msbPos)))
		tMSB[i] = rtypeBit.FromUint64(uint64(tBits[0][i].GetBit(msbPos)))
		MSB[i] = carryMSB[i].Add(tMSB[i])
		if pid == mpcObj.GetHubPid() {
			MSB[i] = MSB[i].Add(aMSB[i])
		}
	}
	// {
	// 	if pid > 0 {
	// 		log.LLvl1("pid", pid, "\nCarryMSB", mpcObj.RevealSymVec(carryMSB)[:5], "\naMSB", aMSB.ToInt()[:5], "\ntMSB",
	// 			mpcObj.RevealSymVec(tMSB)[:5], "\nMSB", mpcObj.RevealSymVec(MSB)[:5])
	// 	}
	// }

	// Flip MSB so that positive result becomes 1
	if pid == mpcObj.GetHubPid() {
		MSB.AddScalar(rtypeBit.FromUint64(1))
	}

	// {
	// 	if pid > 0 {
	// 		log.LLvl1("pid", "\nFlippedMSB", mpcObj.RevealSymVec(MSB)[:5])
	// 	}
	// }

	rConv, rConvBit := mpcObj.ShareRandomBits(rtypeInput, rtypeBit, n, 1, 1, true)
	rConvBit = rConvBit.Transpose()

	// if pid > 0 {
	// 	log.LLvl1("pid", pid, "\nrConv", mpcObj.RevealSymVec(rConv)[:5], "\nrConvBit", mpcObj.RevealSymVec(rConvBit[0])[:5])
	// }

	MSB.Add(rConvBit[0])
	v := mpcObj.RevealSymVec(MSB).ToInt()
	for i := range rConv {
		rConv[i] = rConv[i].Sub(rConv[i].Mul(rtypeInput.FromInt(2 * v[i])))
		if pid == mpcObj.GetHubPid() {
			rConv[i] = rConv[i].Add(rtypeInput.FromInt(v[i]))
		}
	}
	// if pid > 0 {
	// 	log.LLvl1("pid", pid, "\nv", v[:5], "\nOut", mpcObj.RevealSymVec(rConv)[:5])
	// }

	return rConv

	// tBits = mpcObj.RevealSymMat(tBits)
	// sOut := mpc_core.InitRVec(rtype.Zero(), n)
	// for i := 0; i < len(aBits); i++ {
	// 	carry := rtypeBit.Zero()
	// 	sPrint := mpc_core.InitRVec(rtypeBit, numBits)
	// 	for j := numBits - 1; j >= 0; j-- {
	// 		b1 := uint64(aBits[i][j].(mpc_core.LElem2N))
	// 		b2 := uint64(tBits[i][j].(mpc_core.LElem2N))
	// 		c := uint64(carry.(mpc_core.LElem2N))
	// 		sum := b1 + b2 + c
	// 		o1 := rtypeBit.Zero()
	// 		o2 := rtypeBit.Zero()
	// 		if sum == 1 || sum == 3 {
	// 			o1 = rtypeBit.FromInt(1)
	// 		}
	// 		if sum >= 2 {
	// 			o2 = rtypeBit.FromInt(1)
	// 		}
	// 		carry = o2
	// 		if j == 0 && pid == mpcObj.GetHubPid() { // Comparison result
	// 			sOut[i] = rtypeBit.FromInt(1).Sub(o1)
	// 		}
	// 		sPrint[j] = o1
	// 	}
	// 	if i < 5 {
	// 		// xr := mpcObj.RevealSym(x[i])
	// 		// xBits := numToBits(rtypeBit, mpc_core.RVec{xr}, numBits)
	// 		// if pid > 0 {
	// 		// 	log.Print(pid, i, ":", "\ns", sPrint, "\nx", xBits[0], "\na", aBits[i], "\nt", tBits[i])
	// 		// }
	// 	}
	// }
	// return sOut
}

func (mpcObj *MPC) IsPositive(a mpc_core.RVec, binaryVersion bool) mpc_core.RVec {
	if binaryVersion {
		return mpcObj.IsPositive2N(a)
	}

	pid := mpcObj.GetPid()

	n := len(a)
	rtype := a.Type()
	rtypeBit := mpc_core.SElemC(0)
	nbits := rtype.ModBitLength()
	two := rtypeBit.FromInt(2)

	r, rBits := mpcObj.ShareRandomBits(rtype, rtypeBit, n, nbits, 0, false)

	// Compute and reveal 2 * a + r
	c := a.Copy()
	c.MulScalar(rtype.FromInt(2))
	c.Add(r)
	c = mpcObj.RevealSymVec(c)

	cBits := numToBits(rtypeBit, c, nbits)

	noOverflow := mpcObj.LessThanBitsPublic(rBits, cBits)

	cxr := mpc_core.InitRVec(rtypeBit.Zero(), n)
	if pid > 0 {
		for i := range cxr {
			cxr[i] = rBits[i][nbits-1].Sub(two.Mul(cBits[i][nbits-1]).Mul(rBits[i][nbits-1]))
			if pid == 1 {
				cxr[i] = cxr[i].Add(cBits[i][nbits-1])
			}
		}
	}

	lsb := mpcObj.SSMultElemVec(cxr, noOverflow)
	if pid > 0 {
		lsb.MulScalar(two)
		for i := range lsb {
			lsb[i] = lsb[i].Sub(noOverflow[i].Add(cxr[i]))
			if pid == 1 {
				lsb[i] = lsb[i].Add(rtypeBit.One())
			}
		}
	}

	// 0, 1 -> 1, 2
	if pid == 1 {
		lsb.AddScalar(rtypeBit.One())
	}

	return mpcObj.TableLookupWithShareConversion(rtype, lsb, 0)[0]
}

func (mpcObj *MPC) Trunc(a mpc_core.RElem, k, m int) mpc_core.RElem {
	return mpcObj.TruncMat(mpc_core.RMat{mpc_core.RVec{a}}, k, m)[0][0]
}
func (mpcObj *MPC) TruncVec(a mpc_core.RVec, k, m int) mpc_core.RVec {
	return mpcObj.TruncMat(mpc_core.RMat{a}, k, m)[0]
}
func (mpcObj *MPC) TruncMat(a mpc_core.RMat, k, m int) mpc_core.RMat {
	pid := mpcObj.GetPid()
	rtype := a.Type()
	nr, nc := a.Dims()
	out := a.Copy()

	var r, rLow mpc_core.RMat
	if pid == 0 {
		r = mpcObj.Network.Rand.RandMatBits(rtype, nr, nc, rtype.ModBitLength()-2)
		rLow = r.Copy()
		rLow.Trunc(m)

		for p := 1; p < mpcObj.Network.NumParties-1; p++ {
			mpcObj.Network.Rand.SwitchPRG(p)
			rMask := mpcObj.Network.Rand.RandMat(rtype, nr, nc)
			rLowMask := mpcObj.Network.Rand.RandMat(rtype, nr, nc)
			mpcObj.Network.Rand.RestorePRG()
			r.Sub(rMask)
			rLow.Sub(rLowMask)
		}
		mpcObj.Network.SendRData(r, mpcObj.Network.NumParties-1)
		mpcObj.Network.SendRData(rLow, mpcObj.Network.NumParties-1)
	} else if pid == mpcObj.Network.NumParties-1 {
		r = mpcObj.Network.ReceiveRMat(rtype, nr, nc, 0)
		rLow = mpcObj.Network.ReceiveRMat(rtype, nr, nc, 0)
	} else {
		mpcObj.Network.Rand.SwitchPRG(0)
		r = mpcObj.Network.Rand.RandMat(rtype, nr, nc)
		rLow = mpcObj.Network.Rand.RandMat(rtype, nr, nc)
		mpcObj.Network.Rand.RestorePRG()
	}

	// Compute and reveal a + r
	c := a.Copy()
	c.Add(r)
	c = mpcObj.RevealSymMat(c)

	cLow := c.Copy()
	cLow.Trunc(m)

	if pid > 0 {
		out.Add(rLow)
		if pid == 1 {
			out.Sub(cLow)
		}
	}

	key := TypedKey{m, rtype.TypeID()}
	twoInvM, found := mpcObj.invPowCache[key]
	if !found {
		twoInv := rtype.FromInt(2).Inv()
		twoInvM = twoInv
		for p := 0; p < m-1; p++ {
			twoInvM = twoInvM.Mul(twoInv)
		}
		mpcObj.invPowCache[key] = twoInvM
	}

	out.MulScalar(twoInvM)

	return out
}

// Assumes K - F is even
func (mpcObj *MPC) SqrtAndSqrtInverse(a mpc_core.RVec, binaryVersion bool) (b, bInv mpc_core.RVec) {
	if len(a) > mpcObj.divSqrtMaxLen {
		n := len(a)
		c := make(mpc_core.RVec, n)
		cInv := make(mpc_core.RVec, n)
		for remainder, shift := n, 0; remainder > 0; {
			start := shift
			end := shift + mpcObj.divSqrtMaxLen
			if end > n {
				end = n
			}

			fmt.Printf("MPC sqrt/sqrtInv on large vector (%d-%d / %d)\n", start, end, n)

			out, outInv := mpcObj.SqrtAndSqrtInverse(a[start:end], binaryVersion)
			copy(c[start:end], out)
			copy(cInv[start:end], outInv)

			shift = end
			remainder -= end - start
		}
		return c, cInv
	}

	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType().Zero()
	n := len(a)

	nBitsK := mpcObj.GetDataBits()
	nBitsF := mpcObj.GetFracBits()

	numIter := 2 * int(math.Ceil(math.Log2(float64(nBitsK)/3.5)))

	// Initial approximation: 1 / sqrt(aScaled) ~= 2.9581 - 4 * aScaled + 2 * aScaled^2
	var s, sSqrt mpc_core.RVec
	if binaryVersion {
		s, sSqrt = mpcObj.NormalizerEvenExp2N(a, nBitsK)
	} else {
		s, sSqrt = mpcObj.NormalizerEvenExp(a, nBitsK)
	}

	aScaled := mpcObj.SSMultElemVec(a, s)
	aScaled = mpcObj.TruncVec(aScaled, nBitsK, nBitsK-nBitsF)

	aScaledSq := mpcObj.SSMultElemVec(aScaled, aScaled)
	aScaledSq = mpcObj.TruncVec(aScaledSq, nBitsK, nBitsF)

	scaledEst := mpc_core.InitRVec(rtype.Zero(), n)
	if pid > 0 {
		aScaled.MulScalar(rtype.FromInt(4).Neg())
		aScaledSq.MulScalar(rtype.FromInt(2))
		scaledEst.Add(aScaled)
		scaledEst.Add(aScaledSq)

		if pid == mpcObj.GetHubPid() {
			scaledEst.AddScalar(rtype.FromFloat64(2.9581, nBitsF))
		}
	}
	aScaled, aScaledSq = nil, nil

	h := mpcObj.SSMultElemVec(scaledEst, sSqrt)
	h = mpcObj.TruncVec(h, nBitsK/2+nBitsF+2, (nBitsK-nBitsF)/2+1)

	g := h.Copy()
	g.MulScalar(rtype.FromInt(2))

	g = mpcObj.SSMultElemVec(g, a)
	g = mpcObj.TruncVec(g, nBitsK, nBitsF)

	for it := 0; it < numIter; it++ {
		r := mpcObj.SSMultElemVec(h, g)
		r = mpcObj.TruncVec(r, nBitsK, nBitsF)
		r.MulScalar(rtype.One().Neg())

		if pid == 1 {
			r.AddScalar(rtype.FromFloat64(1.5, nBitsF))
		}

		g = mpcObj.SSMultElemVec(g, r)
		h = mpcObj.SSMultElemVec(h, r)

		g = mpcObj.TruncVec(g, nBitsK, nBitsF)
		h = mpcObj.TruncVec(h, nBitsK, nBitsF)
	}

	bInv = h
	bInv.MulScalar(rtype.FromInt(2))
	b = g
	return b, bInv
}

func (mpcObj *MPC) Divide(a, b mpc_core.RVec, binaryVersion bool) mpc_core.RVec {
	log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("MPC Divide called (n = %d) divSqrtMaxLen %d", len(a), mpcObj.divSqrtMaxLen))

	if len(a) > mpcObj.divSqrtMaxLen {
		n := len(a)
		c := make(mpc_core.RVec, n)
		for remainder, shift := n, 0; remainder > 0; {
			start := shift
			end := shift + mpcObj.divSqrtMaxLen
			if end > n {
				end = n
			}

			log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("MPC division on large vector (%d-%d / %d)", start, end, n))

			copy(c[start:end], mpcObj.Divide(a[start:end], b[start:end], binaryVersion))

			shift = end
			remainder -= end - start
		}
		return c
	}

	pid := mpcObj.GetPid()

	if len(a) != len(b) {
		panic("Input vector lengths do not match")
	}

	n := len(a)
	rtype := a.Type()

	nBitsK := mpcObj.dataBits
	nBitsF := mpcObj.fracBits

	numIter := 2*int(math.Ceil(math.Log2(float64(nBitsK)/3.5))) + 1

	// Initial approximation: 1 / bScaled ~= 5.9430 - 10 * bScaled + 5 * bScaled^2
	var s mpc_core.RVec
	if binaryVersion {
		s, _ = mpcObj.NormalizerEvenExp2N(b, nBitsK)
	} else {
		s, _ = mpcObj.NormalizerEvenExp(b, nBitsK)
	}

	bScaled := mpcObj.SSMultElemVec(b, s)
	bScaled = mpcObj.TruncVec(bScaled, nBitsK, nBitsK-nBitsF)

	bScaledSq := mpcObj.SSMultElemVec(bScaled, bScaled)
	bScaledSq = mpcObj.TruncVec(bScaledSq, nBitsK, nBitsF)

	scaledEst := mpc_core.InitRVec(rtype.Zero(), n)
	if pid > 0 {
		bScaled.MulScalar(rtype.FromInt(10).Neg())
		bScaledSq.MulScalar(rtype.FromInt(5))
		scaledEst.Add(bScaled)
		scaledEst.Add(bScaledSq)

		if pid == 1 {
			scaledEst.AddScalar(rtype.FromFloat64(5.9430, nBitsF))
		}
	}
	bScaled, bScaledSq = nil, nil

	w := mpcObj.SSMultElemVec(scaledEst, s)
	w = mpcObj.TruncVec(w, nBitsK+nBitsF+2, nBitsK-nBitsF)

	x := mpcObj.SSMultElemVec(w, b)
	x = mpcObj.TruncVec(x, nBitsK, nBitsF)

	x.MulScalar(rtype.One().Neg())
	if pid == 1 {
		x.AddScalar(rtype.FromFloat64(1, nBitsF))
	}

	y := mpcObj.SSMultElemVec(w, a)
	y = mpcObj.TruncVec(y, nBitsK, nBitsF)

	for it := 0; it < numIter; it++ {

		xr, xm := mpcObj.BeaverPartitionVec(x)
		yr, ym := mpcObj.BeaverPartitionVec(y)

		xpr := xr.Copy()
		if pid > 0 {
			xpr.AddScalar(rtype.FromFloat64(1, nBitsF))
		}

		y = mpcObj.BeaverMultElemVec(yr, ym, xpr, xm)
		x = mpcObj.BeaverMultElemVec(xr, xm, xr, xm)

		y = mpcObj.BeaverReconstructVec(y)
		x = mpcObj.BeaverReconstructVec(x)

		y = mpcObj.TruncVec(y, nBitsK, nBitsF)
		x = mpcObj.TruncVec(x, nBitsK, nBitsF)
	}

	if pid == 1 {
		x.AddScalar(rtype.FromFloat64(1, nBitsF))
	}

	c := mpcObj.SSMultElemVec(x, y)
	c = mpcObj.TruncVec(c, nBitsK, nBitsF)
	return c
}

// func (mpcObj *MPC)LessOrGreaterThan(a mpc_core.RVec, val mpc_core.RElem, isPublic, notLessThan bool) mpc_core.RVec {
// 	valVec := mpc_core.InitRVec(val, len(a))
// 	return LessOrGreaterThanVec(a, valVec, isPublic, notLessThan)
// }

// //LessOrGreaterThan returns vector 0/1 based on if val is less, greater
// func (mpcObj *MPC)LessOrGreaterThanVec(a, val mpc_core.RVec, isPublic, notLessThan bool) mpc_core.RVec {
// 	aCopy := a.Copy()
// 	if pid > 0 {
// 		if isPublic {
// 			if pid == 1 {
// 				aCopy.Sub(val)
// 			}
// 		} else {
// 			aCopy.Sub(val)
// 		}

// 	}

// 	aCopy = IsPositive(aCopy) //0,1
// 	if notLessThan {
// 		//want to know if val is greater than or equal to
// 		return aCopy
// 	}

// 	return mpcObj.FlipBit(aCopy) //call neg on all values?
// }

func (mpcObj *MPC) FlipBit(a mpc_core.RVec) mpc_core.RVec {
	pid := mpcObj.GetPid()

	b := mpc_core.InitRVec(a[0].Zero(), len(a))

	if pid != 0 {
		b.Sub(a)
		if pid == 1 {
			b.AddScalar(a[0].FromInt(1))
		}
	}
	return b

}

func (mpcObj *MPC) NotLessThan(a, b mpc_core.RVec, binaryVersion bool) mpc_core.RVec {
	res := mpcObj.LessThan(a, b, binaryVersion)
	res = mpcObj.FlipBit(res)
	return res
}

func (mpcObj *MPC) NotLessThanPublic(a mpc_core.RVec, bpub mpc_core.RElem, binaryVersion bool) mpc_core.RVec {
	res := mpcObj.LessThanPublic(a, bpub, binaryVersion)
	res = mpcObj.FlipBit(res)
	return res
}

func (mpcObj *MPC) LessThan(a, b mpc_core.RVec, binaryVersion bool) mpc_core.RVec {
	pid := mpcObj.GetPid()

	acopy := a.Copy()
	if pid > 0 {
		acopy.Sub(b)
	}

	res := mpcObj.IsPositive(acopy, binaryVersion)
	res = mpcObj.FlipBit(res)

	return res
}

func (mpcObj *MPC) LessThanPublic(a mpc_core.RVec, bpub mpc_core.RElem, binaryVersion bool) mpc_core.RVec {
	pid := mpcObj.GetPid()

	acopy := a.Copy()
	if pid > 0 {
		if pid == 1 {
			bvec := mpc_core.InitRVec(bpub, len(a))
			acopy.Sub(bvec)
		}
	}

	res := mpcObj.IsPositive(acopy, binaryVersion)

	res = mpcObj.FlipBit(res)
	return res
}

/* Linear algebra routines for eigen-decomposition */
func (mpcObj *MPC) Householder(x mpc_core.RVec) mpc_core.RVec {
	n := len(x)
	rtype := x.Type()
	nBitsK := mpcObj.dataBits
	nBitsF := mpcObj.fracBits
	isBinary := mpcObj.useBooleanShares

	xr, xm := mpcObj.BeaverPartitionVec(x)

	z := mpcObj.BeaverMultElemVec(xr, xm, xr, xm)
	xdot := rtype.Zero()
	for i := range z {
		xdot = xdot.Add(z[i])
	}
	xdot = mpcObj.BeaverReconstruct(xdot)
	xdot = mpcObj.Trunc(xdot, nBitsK, nBitsF)

	xnorm, _ := mpcObj.SqrtAndSqrtInverse(mpc_core.RVec{xdot}, isBinary)

	x1sign := mpcObj.IsPositive(mpc_core.RVec{x[0]}, isBinary)
	x1sign.MulScalar(rtype.FromInt(2))
	if mpcObj.GetPid() == mpcObj.GetHubPid() {
		x1sign.AddScalar(rtype.FromInt(-1))
	}

	shift := mpcObj.SSMultElemVec(xnorm, x1sign)[0]

	sr, sm := mpcObj.BeaverPartition(shift)
	dotShift := mpcObj.BeaverMult(xr[0], xm[0], sr, sm)
	dotShift = mpcObj.BeaverReconstruct(dotShift)
	dotShift = mpcObj.Trunc(dotShift, nBitsK, nBitsF)

	vdot := xdot.Add(dotShift).Mul(rtype.FromInt(2))

	_, vnormInv := mpcObj.SqrtAndSqrtInverse(mpc_core.RVec{vdot}, isBinary)

	invr, invm := mpcObj.BeaverPartition(vnormInv[0])
	vr := xr.Copy()
	vm := xm.Copy()
	vr[0] = vr[0].Add(sr)
	vm[0] = vm[0].Add(sm)

	v := mpc_core.InitRVec(rtype.Zero(), n)
	for i := range v {
		v[i] = mpcObj.BeaverMult(vr[i], vm[i], invr, invm)
	}
	v = mpcObj.BeaverReconstructVec(v)
	return mpcObj.TruncVec(v, nBitsK, nBitsF)
}

func (mpcObj *MPC) QRFactSquare(A mpc_core.RMat) (Q, R mpc_core.RMat) {
	pid := mpcObj.GetPid()

	nBitsK := mpcObj.dataBits
	nBitsF := mpcObj.fracBits

	if len(A) != len(A[0]) {
		panic("QRFactSquare: input matrix is not square")
	}

	n := len(A)
	rtype := A.Type()

	one := rtype.FromFloat64(1, mpcObj.fracBits)

	R = mpc_core.InitRMat(rtype.Zero(), n, n)
	Ap := A.Copy()

	for i := 0; i < n-1; i++ {
		v := make(mpc_core.RMat, 1)
		v[0] = mpcObj.Householder(Ap[0])

		vt := v.Transpose()

		P := mpcObj.SSMultMat(vt, v)
		P = mpcObj.TruncMat(P, nBitsK, nBitsF)

		if pid > 0 {
			P.MulScalar(rtype.FromInt(-2))
			if pid == mpcObj.GetHubPid() {
				for j := range P {
					P[j][j] = P[j][j].Add(one)
				}
			}
		}

		var B mpc_core.RMat
		if i == 0 {
			Q = P.Copy()
			B = mpcObj.SSMultMat(Ap, P)
			B = mpcObj.TruncMat(B, nBitsK, nBitsF)
		} else {
			Qsub := mpc_core.InitRMat(rtype.Zero(), n-i, n)
			if pid > 0 {
				for j := range Qsub {
					Qsub[j] = Q[j+i].Copy()
				}
			}

			// TODO: parallelize
			r0 := mpcObj.SSMultMat(P, Qsub)
			r1 := mpcObj.SSMultMat(Ap, P)
			r0 = mpcObj.TruncMat(r0, nBitsK, nBitsF)
			r1 = mpcObj.TruncMat(r1, nBitsK, nBitsF)

			if pid > 0 {
				for j := range Qsub {
					Q[j+i] = r0[j]
				}
			}
			B = r1
		}

		if pid > 0 {
			for j := 0; j < n-i; j++ {
				R[i+j][i] = B[j][0]
			}
			if i == n-2 {
				R[n-1][n-1] = B[1][1]
			}
		}

		Ap = mpc_core.InitRMat(rtype.Zero(), n-i-1, n-i-1)

		if pid > 0 {
			for j := range Ap {
				for k := range Ap[j] {
					Ap[j][k] = B[j+1][k+1]
				}
			}
		}
	}
	return
}

func (mpcObj *MPC) Tridiag(A mpc_core.RMat) (T, Q mpc_core.RMat) {
	pid := mpcObj.GetPid()

	nBitsK := mpcObj.dataBits
	nBitsF := mpcObj.fracBits

	if len(A) != len(A[0]) {
		panic("Tridiag: input matrix is not square")
	}

	n := len(A)
	rtype := A.Type()

	one := rtype.FromFloat64(1, mpcObj.fracBits)

	T = mpc_core.InitRMat(rtype.Zero(), n, n)
	Q = mpc_core.InitRMat(rtype.Zero(), n, n)
	if pid == mpcObj.GetHubPid() {
		for i := range Q {
			Q[i][i] = one
		}
	}

	Ap := A.Copy()

	for i := 0; i < n-2; i++ {
		x := mpc_core.InitRVec(rtype.Zero(), len(Ap[0])-1)
		if pid > 0 {
			for j := range x {
				x[j] = Ap[0][j+1]
			}
		}

		v := make(mpc_core.RMat, 1)
		v[0] = mpcObj.Householder(x)
		vt := v.Transpose()

		vv := mpcObj.SSMultMat(vt, v)
		vv = mpcObj.TruncMat(vv, nBitsK, nBitsF)

		P := mpc_core.InitRMat(rtype.Zero(), len(Ap), len(Ap))
		if pid > 0 {
			if pid == mpcObj.GetHubPid() {
				P[0][0] = one
			}
			for j := 1; j < len(Ap); j++ {
				for k := 1; k < len(Ap); k++ {
					P[j][k] = vv[j-1][k-1].Mul(rtype.FromInt(-2))
					if pid == mpcObj.GetHubPid() && j == k {
						P[j][k] = P[j][k].Add(one)
					}
				}
			}
		}

		PAp := mpcObj.SSMultMat(P, Ap)
		PAp = mpcObj.TruncMat(PAp, nBitsK, nBitsF)

		B := mpcObj.SSMultMat(PAp, P)
		B = mpcObj.TruncMat(B, nBitsK, nBitsF)

		Qsub := mpc_core.InitRMat(rtype.Zero(), n, n-i)
		if pid > 0 {
			for j := range Qsub {
				for k := range Qsub[j] {
					Qsub[j][k] = Q[j][k+i]
				}
			}
		}

		Qsub = mpcObj.SSMultMat(Qsub, P)
		Qsub = mpcObj.TruncMat(Qsub, nBitsK, nBitsF)
		if pid > 0 {
			for j := range Qsub {
				for k := range Qsub[j] {
					Q[j][k+i] = Qsub[j][k]
				}
			}
		}

		if pid > 0 {
			T[i][i] = B[0][0]
			T[i+1][i] = B[1][0]
			T[i][i+1] = B[0][1]
			if i == n-3 {
				T[i+1][i+1] = B[1][1]
				T[i+1][i+2] = B[1][2]
				T[i+2][i+1] = B[2][1]
				T[i+2][i+2] = B[2][2]
			}
		}

		Ap = mpc_core.InitRMat(rtype.Zero(), len(B)-1, len(B[0])-1)
		if pid > 0 {
			for j := range Ap {
				for k := range Ap[j] {
					Ap[j][k] = B[j+1][k+1]
				}
			}
		}
	}
	return
}

func (mpcObj *MPC) Swap(v1 mpc_core.RVec, v2 mpc_core.RVec, isFlip mpc_core.RElem) (o1 mpc_core.RVec, o2 mpc_core.RVec) {
	// o1 = (1 - isFlip) * v1 + isFlip * v2
	//    = v1 + isFlip * (v2 - v1)
	// o2 = v1 + v2 - o1

	d := v2.Copy()
	d.Sub(v1)

	m := mpcObj.SSMultElemVecScalar(d, isFlip)

	o1 = v1.Copy()
	o1.Add(m)

	o2 = v2.Copy()
	o2.Sub(m)

	return
}

func (mpcObj *MPC) SortRowsDescend(A mpc_core.RMat, w mpc_core.RVec) (ASorted mpc_core.RMat, wSorted mpc_core.RVec) {
	ASorted = A.Copy()
	W := mpc_core.RMat{w}.Transpose()

	n := len(ASorted)
	isBinary := mpcObj.useBooleanShares

	for i := 0; i < n-1; i++ {
		for j := n - 1; j > i; j-- {
			// Compare W[j] and W[j-1] and swap if W[j] > W[j-1]
			isFlip := mpcObj.IsPositive(mpc_core.RVec{W[j][0].Sub(W[j-1][0])}, isBinary)[0]
			ASorted[j], ASorted[j-1] = mpcObj.Swap(ASorted[j], ASorted[j-1], isFlip)
			W[j], W[j-1] = mpcObj.Swap(W[j], W[j-1], isFlip)
		}
	}

	wSorted = W.Transpose()[0]
	return
}

func (mpcObj *MPC) EigenDecomp(A mpc_core.RMat) (V mpc_core.RMat, L mpc_core.RVec) {
	pid := mpcObj.GetPid()

	const ITER_PER_EVAL int = 5
	nBitsK := mpcObj.dataBits
	nBitsF := mpcObj.fracBits

	if len(A) != len(A[0]) {
		panic("EigenDecomp: input matrix is not square")
	}

	n := len(A)
	rtype := A.Type()
	L = mpc_core.InitRVec(rtype.Zero(), n)

	Ap, Q := mpcObj.Tridiag(A)
	if pid == 0 {
		V = mpc_core.InitRMat(rtype.Zero(), n, n)
	} else {
		V = Q.Transpose()
	}

	for i := n - 1; i >= 1; i-- {
		fmt.Println(fmt.Sprintf("EigenDecomp: %d-th eigenvalue", i))

		for it := 0; it < ITER_PER_EVAL; it++ {
			shift := Ap[i][i]
			if pid > 0 {
				for j := range Ap {
					Ap[j][j] = Ap[j][j].Sub(shift)
				}
			}

			Q, R := mpcObj.QRFactSquare(Ap)
			Ap = mpcObj.SSMultMat(Q, R)
			Ap = mpcObj.TruncMat(Ap, nBitsK, nBitsF)

			if pid > 0 {
				for j := range Ap {
					Ap[j][j] = Ap[j][j].Add(shift)
				}
			}

			Vsub := mpc_core.InitRMat(rtype.Zero(), i+1, n)
			if pid > 0 {
				for j := range Vsub {
					for k := range Vsub[j] {
						Vsub[j][k] = V[j][k]
					}
				}
			}

			Vsub = mpcObj.SSMultMat(Q, Vsub)
			Vsub = mpcObj.TruncMat(Vsub, nBitsK, nBitsF)

			if pid > 0 {
				for j := range Vsub {
					for k := range Vsub[j] {
						V[j][k] = Vsub[j][k]
					}
				}
			}
		}

		L[i] = Ap[i][i]
		if i == 1 {
			L[0] = Ap[0][0]
		}

		newAp := mpc_core.InitRMat(rtype.Zero(), i, i)
		if pid > 0 {
			for j := range newAp {
				for k := range newAp[j] {
					newAp[j][k] = Ap[j][k]
				}
			}
		}
		Ap = newAp
	}

	fmt.Println("EigenDecomp: complete")
	return
}

/* PARALLEL ROUTINES*/
func (mpcObjs ParallelMPC) DisableLogging() {
	for i := range mpcObjs {
		mpcObjs[i].Network.DisableLogging()
	}
}

func (mpcObjs ParallelMPC) EnableLogging() {
	for i := range mpcObjs {
		mpcObjs[i].Network.EnableLogging()
	}
}

func (mpcObjs ParallelMPC) GetNetworks() ParallelNetworks {
	net := make(ParallelNetworks, len(mpcObjs))
	for i := range net {
		net[i] = mpcObjs[i].Network
	}
	return net
}

type MpcRoutine func(*MPC, mpc_core.RMat, mpc_core.RElem) mpc_core.RVec

func (mpcObjs ParallelMPC) NotLessThanPublic(a mpc_core.RVec, b mpc_core.RElem, binaryVersion bool) mpc_core.RVec {
	return mpcObjs.runParallel(mpc_core.RMat{a}, b, "NotLessThanPublic", mpcObjs[0].divSqrtMaxLen,
		func(mpc *MPC, mat mpc_core.RMat, aux mpc_core.RElem) mpc_core.RVec {
			return mpc.NotLessThanPublic(mat[0], aux, binaryVersion)
		})
}

func (mpcObjs ParallelMPC) RevealSymVec(a mpc_core.RVec) mpc_core.RVec {
	return mpcObjs.runParallel(mpc_core.RMat{a}, nil, "RevealSymVec", mpcObjs[0].divSqrtMaxLen,
		func(mpc *MPC, mat mpc_core.RMat, aux mpc_core.RElem) mpc_core.RVec {
			return mpc.RevealSymVec(mat[0])
		})
}

func (mpcObjs ParallelMPC) SSSquareElemVec(a mpc_core.RVec) mpc_core.RVec {
	return mpcObjs.runParallel(mpc_core.RMat{a}, nil, "SSSquareElemVec", mpcObjs[0].divSqrtMaxLen,
		func(mpc *MPC, mat mpc_core.RMat, aux mpc_core.RElem) mpc_core.RVec {
			return mpc.SSSquareElemVec(mat[0])
		})
}

func (mpcObjs ParallelMPC) SSMultElemVec(a, b mpc_core.RVec) mpc_core.RVec {
	return mpcObjs.runParallel(mpc_core.RMat{a, b}, nil, "SSMultElemVec", mpcObjs[0].divSqrtMaxLen,
		func(mpc *MPC, mat mpc_core.RMat, aux mpc_core.RElem) mpc_core.RVec {
			return mpc.SSMultElemVec(mat[0], mat[1])
		})
}

func (mpcObjs ParallelMPC) LessThan(a, b mpc_core.RVec, binaryVersion bool) mpc_core.RVec {
	return mpcObjs.runParallel(mpc_core.RMat{a, b}, nil, "LessThan", mpcObjs[0].divSqrtMaxLen,
		func(mpc *MPC, mat mpc_core.RMat, aux mpc_core.RElem) mpc_core.RVec {
			return mpc.LessThan(mat[0], mat[1], binaryVersion)
		})
}

func (mpcObjs ParallelMPC) Divide(a, b mpc_core.RVec, binaryVersion bool) mpc_core.RVec {
	return mpcObjs.runParallel(mpc_core.RMat{a, b}, nil, "Divide", mpcObjs[0].divSqrtMaxLen,
		func(mpc *MPC, mat mpc_core.RMat, aux mpc_core.RElem) mpc_core.RVec {
			return mpc.Divide(mat[0], mat[1], binaryVersion)
		})
}

func (mpcObjs ParallelMPC) Sqrt(a mpc_core.RVec, binaryVersion bool) mpc_core.RVec {
	return mpcObjs.runParallel(mpc_core.RMat{a}, nil, "Sqrt", mpcObjs[0].divSqrtMaxLen,
		func(mpc *MPC, mat mpc_core.RMat, aux mpc_core.RElem) mpc_core.RVec {
			res, _ := mpc.SqrtAndSqrtInverse(mat[0], binaryVersion)
			return res
		})
}

func (mpcObjs ParallelMPC) SqrtInv(a mpc_core.RVec, binaryVersion bool) mpc_core.RVec {
	return mpcObjs.runParallel(mpc_core.RMat{a}, nil, "SqrtInv", mpcObjs[0].divSqrtMaxLen,
		func(mpc *MPC, mat mpc_core.RMat, aux mpc_core.RElem) mpc_core.RVec {
			_, res := mpc.SqrtAndSqrtInverse(mat[0], binaryVersion)
			return res
		})
}

func (mpcObjs ParallelMPC) IsPositive(a mpc_core.RVec, binaryVersion bool) mpc_core.RVec {
	return mpcObjs.runParallel(mpc_core.RMat{a}, nil, "IsPositive", mpcObjs[0].divSqrtMaxLen,
		func(mpc *MPC, mat mpc_core.RMat, aux mpc_core.RElem) mpc_core.RVec {
			return mpc.IsPositive(mat[0], binaryVersion)
		})
}

// If batchSize > 0 and smaller than len(a[0]), then break up the matrix into chucks to process sequentially (each chunk is parallelized)
func (mpcObjs ParallelMPC) runParallel(a mpc_core.RMat, aux mpc_core.RElem, name string, batchSize int, fn MpcRoutine) mpc_core.RVec {
	n := len(a[0])

	log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("runParallel called (%s): initializing n %d batchSize %d", name, n, batchSize))

	res := mpc_core.InitRVec(mpcObjs[0].GetRType().Zero(), n)

	if batchSize > 0 && n > batchSize {

		for startIndex, endIndex := 0, batchSize; startIndex < n; startIndex, endIndex = startIndex+batchSize, endIndex+batchSize {
			if endIndex > n {
				endIndex = n
			}

			log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("runParallel (%s): working on %d-%d / %d", name, startIndex, endIndex, n))

			aSub := make(mpc_core.RMat, len(a))
			for row := range aSub {
				aSub[row] = a[row][startIndex:endIndex]
			}

			out := mpcObjs.runParallel(aSub, aux, name, 0, fn)

			copy(res[startIndex:endIndex], out)
		}

	} else {

		numThreads := len(mpcObjs)
		batchFloat := float64(n) / float64(numThreads)
		divSqrtMaxLen := mpcObjs[0].divSqrtMaxLen

		var wg sync.WaitGroup
		startIndex, endIndex := 0, 0
		for i := 0; i < numThreads; i++ {
			if i == numThreads-1 {
				endIndex = n
			} else {
				endIndex = int(batchFloat * float64(i+1))
			}

			if endIndex > startIndex {
				aSub := make(mpc_core.RMat, len(a))
				for row := range aSub {
					aSub[row] = a[row][startIndex:endIndex]
				}

				wg.Add(1)
				go func(threadID, startIndex, endIndex, divSqrtMaxLen int, aSub mpc_core.RMat) {
					defer wg.Done()
					mpcObjs[threadID].divSqrtMaxLen = divSqrtMaxLen
					tmp := fn(mpcObjs[threadID], aSub, aux)
					copy(res[startIndex:endIndex], tmp)
					// log.LLvl1(fmt.Sprintf("runParallel (%s): processed %d-%d / %d", name, startIndex, endIndex, n))
				}(i, startIndex, endIndex, divSqrtMaxLen, aSub)
			}

			startIndex = endIndex
		}
		wg.Wait()

	}

	return res
}

// If the column norms of A are too large or small there can be precision issues
// Closer they are to unit norm the better
func (mpcObj *MPC) MatrixInverseSVD(A mpc_core.RMat) (AInv mpc_core.RMat) {

	Ar, Am := mpcObj.BeaverPartitionMat(A)
	AAt := mpcObj.BeaverMultMat(Ar, Am, Ar.Transpose(), Am.Transpose())
	AAt = mpcObj.BeaverReconstructMat(AAt)
	AAt = mpcObj.TruncMat(AAt, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	Ut, S2 := mpcObj.EigenDecomp(AAt)

	_, SInv := mpcObj.SqrtAndSqrtInverse(S2, false)
	SInvr, SInvm := mpcObj.BeaverPartitionVec(SInv)

	Ur, Um := mpcObj.BeaverPartitionMat(Ut.Transpose())
	UtA := mpcObj.BeaverMultMat(Ur.Transpose(), Um.Transpose(), Ar, Am)
	UtA = mpcObj.BeaverReconstructMat(UtA)
	UtA = mpcObj.TruncMat(UtA, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	UtAr, UtAm := mpcObj.BeaverPartitionMat(UtA)

	V := make(mpc_core.RMat, len(UtA))
	for i := range V {
		V[i] = mpcObj.BeaverMultMat(mpc_core.RMat{mpc_core.RVec{SInvr[i]}}, mpc_core.RMat{mpc_core.RVec{SInvm[i]}}, mpc_core.RMat{UtAr[i]}, mpc_core.RMat{UtAm[i]})[0]
	}
	V = mpcObj.BeaverReconstructMat(V)
	V = mpcObj.TruncMat(V, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	Vr, Vm := mpcObj.BeaverPartitionMat(V)

	SInvV := make(mpc_core.RMat, len(V))
	for i := range SInvV {
		SInvV[i] = mpcObj.BeaverMultMat(mpc_core.RMat{mpc_core.RVec{SInvr[i]}}, mpc_core.RMat{mpc_core.RVec{SInvm[i]}}, mpc_core.RMat{Vr[i]}, mpc_core.RMat{Vm[i]})[0]
	}
	SInvV = mpcObj.BeaverReconstructMat(SInvV)
	SInvV = mpcObj.TruncMat(SInvV, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	SInvVr, SInvVm := mpcObj.BeaverPartitionMat(SInvV)
	USInvV := mpcObj.BeaverMultMat(Ur, Um, SInvVr, SInvVm)
	USInvV = mpcObj.BeaverReconstructMat(USInvV)
	USInvV = mpcObj.TruncMat(USInvV, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	AInv = USInvV.Transpose()
	return
}

// Assumes A is a symmetric, positive definite matrix (positive eigenvalues)
// Outputs both AInv and AInvSqrt, where AInv = AInvSqrt' * AInvSqrt
func (mpcObj *MPC) MatrixInverseSymPos(ASymPosDef mpc_core.RMat) (AInv, AInvSqrt mpc_core.RMat) {
	Vt, L := mpcObj.EigenDecomp(ASymPosDef)

	// Compute AInvSqrt = L^-(1/2) * Vt
	_, LInvSqrt := mpcObj.SqrtAndSqrtInverse(L, false)

	Vtr, Vtm := mpcObj.BeaverPartitionMat(Vt)
	Lr, Lm := mpcObj.BeaverPartitionVec(LInvSqrt)

	AInvSqrt = make(mpc_core.RMat, len(Vt))
	for i := range Vt {
		AInvSqrt[i] = mpcObj.BeaverMultMat(mpc_core.RMat{mpc_core.RVec{Lr[i]}}, mpc_core.RMat{mpc_core.RVec{Lm[i]}}, mpc_core.RMat{Vtr[i]}, mpc_core.RMat{Vtm[i]})[0]
	}
	AInvSqrt = mpcObj.BeaverReconstructMat(AInvSqrt)

	AInvSqrt = mpcObj.TruncMat(AInvSqrt, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	AInv = mpcObj.SSMultMat(AInvSqrt.Transpose(), AInvSqrt)
	AInv = mpcObj.TruncMat(AInv, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	return
}
