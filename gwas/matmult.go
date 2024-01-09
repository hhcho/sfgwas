package gwas

import (
	"fmt"
	"math/bits"
	"runtime"
	"sync"
	"time"
	"unsafe"

	"go.dedis.ch/onet/v3/log"

	"github.com/hhcho/sfgwas-private/crypto"
	"github.com/ldsec/lattigo/v2/ring"

	"math"

	"github.com/hhcho/sfgwas-private/mpc"
	"github.com/ldsec/lattigo/v2/ckks"

	"gonum.org/v1/gonum/mat"
)

// Multiply Q (kp by nsnp) with X (nsnp by nind) with lazy normalization of X
// Compute Q * S * (X - m * 1^T) as (Q * S) * X - ((Q * S) * m) * 1^T
// S: diagonal matrix containing 1/stdev of each SNP, m: column vector containing mean of each SNP
func QXLazyNormStream(cps *crypto.CryptoParams, mpcObj *mpc.MPC, Q crypto.CipherMatrix, Xcachefile string, XMean, XStdInv crypto.CipherVector, numInd int) (out crypto.CipherMatrix) {
	if mpcObj.GetPid() == 0 {
		return
	}
	slots := cps.GetSlots()

	start := time.Now()

	// Compute Q * S
	QS := make(crypto.CipherMatrix, len(Q))
	for i := range Q {
		QS[i] = crypto.CMult(cps, Q[i], XStdInv)
	}

	// Compute (Q * S) * X
	out = MatMult4StreamCompute(cps, QS, 5, Xcachefile)

	out = mpcObj.Network.BootstrapMatAll(cps, out) // TODO

	// Compute (Q * S) * m
	QSm := make(crypto.CipherVector, len(Q))
	for i := range Q {
		QSm[i] = crypto.InnerProd(cps, QS[i], XMean) // Already has output value in all slots
	}

	// Compute (Q * S) * X - ((Q * S) * m) * 1^T
	for i := range QS {
		cps.WithEvaluator(func(eval ckks.Evaluator) error {
			for j := range out[i] {
				eval.Sub(out[i][j], QSm[i], out[i][j])
			}
			return nil
		})

		// Zero out empty slots at the end
		// TODO: Is there a better way?
		for j := range out[i] {
			var N int
			if j < len(out[i])-1 {
				N = slots
			} else {
				N = ((numInd - 1) % slots) + 1
			}
			out[i][j] = crypto.MaskTrunc(cps, out[i][j], N)
		}
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Matrix multiplication complete,", time.Since(start))

	return
}

// Multiply Q (kp by nind) with X^T (nind by nsnp) with lazy normalization of X
// Compute Q * (X^T - 1 * m^T) * S as ((Q * X^T) - ((Q * 1) * m^T)) * S
// S: diagonal matrix containing 1/stdev of each SNP, m: column vector containing mean of each SNP
// TODO: multiply with S AFTER aggregation across parties, that way bootstrap once for all
func QXtLazyNormStream(cps *crypto.CryptoParams, mpcObj *mpc.MPC, Q crypto.CipherMatrix, XTcachefile string, XMean, XStdInv crypto.CipherVector) (out crypto.CipherMatrix) {
	if mpcObj.GetPid() == 0 {
		return
	}

	start := time.Now()

	// Compute Q * X^T
	out = MatMult4StreamCompute(cps, Q, 5, XTcachefile)
	out = mpcObj.Network.BootstrapMatAll(cps, out) // TODO change to bootstrapping after aggregation

	// Compute (Q * X^T) - ((Q * 1) * m^T)
	for i := range out {
		rowSum := crypto.InnerSumAll(cps, Q[i])

		Q1m := crypto.CMultScalar(cps, XMean, rowSum)

		cps.WithEvaluator(func(eval ckks.Evaluator) error {
			for j := range out[i] {
				eval.Sub(out[i][j], Q1m[j], out[i][j])
			}
			return nil
		})
	}

	// Compute ((Q * X^T) - ((Q * 1) * m^T)) * S
	for i := range out {
		out[i] = crypto.CMult(cps, out[i], XStdInv)
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Matrix multiplication complete,", time.Since(start))

	return
}

type matmulPlainInnerFn func(*crypto.CryptoParams, crypto.CipherVector, crypto.PlainMatrix, int) crypto.CipherVector
type matmulInnerFn func(*crypto.CryptoParams, crypto.CipherVector, crypto.CipherMatrix, int) crypto.CipherVector

func DCMatMulAAtB(cryptoParams *crypto.CryptoParams, mpcObj *mpc.MPC, A crypto.CipherMatrix, B crypto.CipherMatrix,
	nrows []int, ncol_out int, innerFn matmulInnerFn) crypto.CipherMatrix {

	slots := cryptoParams.GetSlots()
	pid := mpcObj.GetPid()

	out := crypto.CZeroMat(cryptoParams, ((nrows[pid]-1)/slots)+1, ncol_out)

	cTQloc := make(crypto.CipherVector, ncol_out)
	for c := range A {
		var wg sync.WaitGroup
		for j := range cTQloc {
			wg.Add(1)
			go func(j int) {
				defer wg.Done()
				innerProd := innerFn(cryptoParams, A[c], B, j)
				cTQloc[j] = crypto.InnerSumAll(cryptoParams, innerProd)
			}(j)
		}
		wg.Wait()

		cTQ := mpcObj.Network.AggregateCVec(cryptoParams, cTQloc)

		for j := 0; j < ncol_out; j++ {
			wg.Add(1)
			go func(j int) {
				defer wg.Done()
				ccTQ := crypto.CMult(cryptoParams, A[c], crypto.CipherVector{cTQ[j]})
				out[j] = crypto.CAdd(cryptoParams, out[j], ccTQ)
			}(j)
		}
		wg.Wait()
	}

	return out
}

func DCMatMulAAtBPlain(cryptoParams *crypto.CryptoParams, mpcObj *mpc.MPC, A crypto.CipherMatrix, B crypto.PlainMatrix,
	nrows []int, ncol_out int, innerFn matmulPlainInnerFn) crypto.CipherMatrix {

	slots := cryptoParams.GetSlots()
	pid := mpcObj.GetPid()

	out := crypto.CZeroMat(cryptoParams, int((nrows[pid]-1)/slots)+1, ncol_out)

	var wg sync.WaitGroup
	for c := range A {
		cTQloc := make(crypto.CipherVector, ncol_out)

		for j := range cTQloc {
			wg.Add(1)
			go func(j int) {
				defer wg.Done()
				innerProd := innerFn(cryptoParams, A[c], B, j)
				cTQloc[j] = crypto.InnerSumAll(cryptoParams, innerProd)
			}(j)
		}
		wg.Wait()

		cTQ := mpcObj.Network.AggregateCVec(cryptoParams, cTQloc)

		for j := 0; j < ncol_out; j++ {
			wg.Add(1)
			go func(j int) {
				defer wg.Done()
				ccTQ := crypto.CMult(cryptoParams, A[c], crypto.CipherVector{cTQ[j]})
				out[j] = crypto.CAdd(cryptoParams, out[j], ccTQ)
			}(j)
		}
		wg.Wait()
	}

	return out
}

type uint128 struct {
	hi uint64
	lo uint64
}

type CipherAccV1 struct {
	acc00 []uint128
	acc01 []uint128
	acc10 []uint128
	acc11 []uint128
}

type CipherAccV2 struct {
	acc0 [][]uint128
	acc1 [][]uint128
}

type CipherVectorAccV1 []CipherAccV1

type CipherVectorAccV2 struct {
	val []CipherAccV2
}

func NewCipherVectorAccV1(cryptoParams *crypto.CryptoParams, n int) CipherVectorAccV1 {
	N := cryptoParams.Params.N()
	out := make(CipherVectorAccV1, n)
	for i := range out {
		out[i].acc00 = make([]uint128, N)
		out[i].acc01 = make([]uint128, N)
		out[i].acc10 = make([]uint128, N)
		out[i].acc11 = make([]uint128, N)
	}
	return out
}

func NewCipherVectorAccV2(cryptoParams *crypto.CryptoParams, n int, level int) CipherVectorAccV2 {
	N := cryptoParams.Params.N()
	var out CipherVectorAccV2
	out.val = make([]CipherAccV2, n)
	for i := range out.val {
		out.val[i].acc0 = make([][]uint128, level)
		out.val[i].acc1 = make([][]uint128, level)
		for l := 0; l < level; l++ {
			out.val[i].acc0[l] = make([]uint128, N)
			out.val[i].acc1[l] = make([]uint128, N)
		}
	}

	return out
}

func MulCoeffsAndAdd128(a, b []uint64, c []uint128) {

	var hi, lo, carry uint64

	for j := 0; j < len(a); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&a[j]))
		y := (*[8]uint64)(unsafe.Pointer(&b[j]))
		z := (*[8]uint128)(unsafe.Pointer(&c[j]))

		hi, lo = bits.Mul64(x[0], y[0])
		z[0].lo, carry = bits.Add64(z[0].lo, lo, 0)
		z[0].hi += hi + carry

		hi, lo = bits.Mul64(x[1], y[1])
		z[1].lo, carry = bits.Add64(z[1].lo, lo, 0)
		z[1].hi += hi + carry

		hi, lo = bits.Mul64(x[2], y[2])
		z[2].lo, carry = bits.Add64(z[2].lo, lo, 0)
		z[2].hi += hi + carry

		hi, lo = bits.Mul64(x[3], y[3])
		z[3].lo, carry = bits.Add64(z[3].lo, lo, 0)
		z[3].hi += hi + carry

		hi, lo = bits.Mul64(x[4], y[4])
		z[4].lo, carry = bits.Add64(z[4].lo, lo, 0)
		z[4].hi += hi + carry

		hi, lo = bits.Mul64(x[5], y[5])
		z[5].lo, carry = bits.Add64(z[5].lo, lo, 0)
		z[5].hi += hi + carry

		hi, lo = bits.Mul64(x[6], y[6])
		z[6].lo, carry = bits.Add64(z[6].lo, lo, 0)
		z[6].hi += hi + carry

		hi, lo = bits.Mul64(x[7], y[7])
		z[7].lo, carry = bits.Add64(z[7].lo, lo, 0)
		z[7].hi += hi + carry
	}
}

func ReduceAndAddUint128(in []uint128, out []uint64, qInv, q uint64) {

	var hhi uint64

	for j := 0; j < len(in); j = j + 8 {

		x := (*[8]uint128)(unsafe.Pointer(&in[j]))
		y := (*[8]uint64)(unsafe.Pointer(&out[j]))

		hhi, _ = bits.Mul64(x[0].lo*qInv, q)
		y[0] += x[0].hi - hhi + q

		hhi, _ = bits.Mul64(x[1].lo*qInv, q)
		y[1] += x[1].hi - hhi + q

		hhi, _ = bits.Mul64(x[2].lo*qInv, q)
		y[2] += x[2].hi - hhi + q

		hhi, _ = bits.Mul64(x[3].lo*qInv, q)
		y[3] += x[3].hi - hhi + q

		hhi, _ = bits.Mul64(x[4].lo*qInv, q)
		y[4] += x[4].hi - hhi + q

		hhi, _ = bits.Mul64(x[5].lo*qInv, q)
		y[5] += x[5].hi - hhi + q

		hhi, _ = bits.Mul64(x[6].lo*qInv, q)
		y[6] += x[6].hi - hhi + q

		hhi, _ = bits.Mul64(x[7].lo*qInv, q)
		y[7] += x[7].hi - hhi + q
	}
}

func ModularReduceV1(cryptoParams *crypto.CryptoParams, cva CipherVectorAccV1, outScale float64) crypto.CipherVector {
	N := cryptoParams.Params.N()
	ringQ, _ := ring.NewRing(N, cryptoParams.Params.Qi())

	out := make(crypto.CipherVector, len(cva))
	for i := range out {
		ct := ckks.NewCiphertext(cryptoParams.Params, 1, 1, outScale)
		ReduceAndAddUint128(cva[i].acc00, ct.Value()[0].Coeffs[0], ringQ.MredParams[0], ringQ.Modulus[0])
		ReduceAndAddUint128(cva[i].acc01, ct.Value()[1].Coeffs[0], ringQ.MredParams[0], ringQ.Modulus[0])
		ReduceAndAddUint128(cva[i].acc10, ct.Value()[0].Coeffs[1], ringQ.MredParams[1], ringQ.Modulus[1])
		ReduceAndAddUint128(cva[i].acc11, ct.Value()[1].Coeffs[1], ringQ.MredParams[1], ringQ.Modulus[1])
		out[i] = ct
	}

	return out
}

func ModularReduceV2(cryptoParams *crypto.CryptoParams, cva CipherVectorAccV2, outScale float64) crypto.CipherVector {
	N := cryptoParams.Params.N()
	ringQ, _ := ring.NewRing(N, cryptoParams.Params.Qi())
	level := len(cva.val[0].acc0)

	out := make(crypto.CipherVector, len(cva.val))
	for i := range out {
		ct := ckks.NewCiphertext(cryptoParams.Params, 1, level-1, outScale)
		for l := 0; l < level; l++ {
			mredParams := ringQ.MredParams[l]
			qi := ringQ.Modulus[l]
			ReduceAndAddUint128(cva.val[i].acc0[l], ct.Value()[0].Coeffs[l], mredParams, qi)
			ReduceAndAddUint128(cva.val[i].acc1[l], ct.Value()[1].Coeffs[l], mredParams, qi)
		}
		err := cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			return eval.Reduce(ct, ct)
		})
		if err != nil {
			panic(err)
		}
		out[i] = ct
	}
	return out
}

// Multiply X and Y to add to Acc without modular reduction
func CPMultAccWithoutMRedV1(cryptoParams *crypto.CryptoParams, X crypto.CipherVector, Y crypto.PlainVector, Acc CipherVectorAccV1) {
	for i := range X {
		if X[i] != nil && Y[i] != nil {
			MulCoeffsAndAdd128(X[i].Value()[0].Coeffs[0], Y[i].Value()[0].Coeffs[0], Acc[i].acc00)
			MulCoeffsAndAdd128(X[i].Value()[1].Coeffs[0], Y[i].Value()[0].Coeffs[0], Acc[i].acc01)
			MulCoeffsAndAdd128(X[i].Value()[0].Coeffs[1], Y[i].Value()[0].Coeffs[1], Acc[i].acc10)
			MulCoeffsAndAdd128(X[i].Value()[1].Coeffs[1], Y[i].Value()[0].Coeffs[1], Acc[i].acc11)
		}
	}
}

func CPMultAccWithoutMRedV2(X crypto.CipherVector, Y crypto.PlainVector, Acc CipherVectorAccV2) {
	n := len(Acc.val)
	for i := 0; i < n; i++ {
		// Broadcasting
		xi, yi := i, i
		if len(X) == 1 {
			xi = 0
		}
		if len(Y) == 1 {
			yi = 0
		}

		if X[xi] != nil && Y[yi] != nil {
			for l := 0; l < len(Acc.val[i].acc0); l++ {
				MulCoeffsAndAdd128(X[xi].Value()[0].Coeffs[l], Y[yi].Value()[0].Coeffs[l], Acc.val[i].acc0[l])
				MulCoeffsAndAdd128(X[xi].Value()[1].Coeffs[l], Y[yi].Value()[0].Coeffs[l], Acc.val[i].acc1[l])
			}
		}
	}
}

func ToMontgomeryForm(cryptoParams *crypto.CryptoParams, pt crypto.PlainVector) {
	N := cryptoParams.Params.N()
	ringQ, _ := ring.NewRing(N, cryptoParams.Params.Qi())
	for i := range pt {
		if pt[i] != nil {
			MFormLvl(ringQ, pt[i].Level(), pt[i].Value()[0], pt[i].Value()[0])
		}
	}
}

func MFormLvl(r *ring.Ring, level int, p1, p2 *ring.Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		bredParams := r.BredParams[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MForm(x[0], qi, bredParams)
			z[1] = MForm(x[1], qi, bredParams)
			z[2] = MForm(x[2], qi, bredParams)
			z[3] = MForm(x[3], qi, bredParams)
			z[4] = MForm(x[4], qi, bredParams)
			z[5] = MForm(x[5], qi, bredParams)
			z[6] = MForm(x[6], qi, bredParams)
			z[7] = MForm(x[7], qi, bredParams)
		}
	}
}

func MForm(a, q uint64, u []uint64) (r uint64) {
	mhi, _ := bits.Mul64(a, u[1])
	r = -(a*u[0] + mhi) * q
	if r >= q {
		r -= q
	}
	return
}

// Cheat to compute ground truth
func CPMatMult0(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, B *mat.Dense) crypto.CipherMatrix {
	n, m := B.Dims()
	s := len(A)
	Araw := crypto.DecryptFloatMatrix(cryptoParams, A, n)
	out := make([][]float64, s)
	for i := range Araw {
		vec := mat.NewVecDense(n, Araw[i])
		vec.MulVec(B.T(), vec)

		out[i] = make([]float64, m)
		for j := range out[i] {
			out[i][j] = vec.AtVec(j)
		}
	}
	enc, _, _, _ := crypto.EncryptFloatMatrixRow(cryptoParams, out)
	return enc
}

// All-pairwise inner product version
func CPMatMult1(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, B *mat.Dense) crypto.CipherMatrix {
	n, m := B.Dims()
	s := len(A)
	slots := cryptoParams.GetSlots()
	m_ct := int((m-1)/slots) + 1

	out := crypto.CZeroMat(cryptoParams, m_ct, s)

	tmp := make([]float64, n)
	for j := 0; j < m; j++ {
		mat.Col(tmp, j, B)
		pv, _ := crypto.EncodeFloatVector(cryptoParams, tmp)

		ctid := int(j / slots)
		slotid := j % slots

		for i := 0; i < s; i++ {
			prod := crypto.CPMult(cryptoParams, A[i], pv)
			ct := crypto.InnerSumAll(cryptoParams, prod)
			ct = crypto.Mask(cryptoParams, ct, slotid, false)
			out[i][ctid] = crypto.Add(cryptoParams, out[i][ctid], ct)
		}
	}
	return out
}

// Duplicate individual elements in A
func CPMatMult2(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, B *mat.Dense) crypto.CipherMatrix {
	n, m := B.Dims()
	s := len(A)
	slots := cryptoParams.GetSlots()
	m_ct := int((m-1)/slots) + 1

	out := crypto.CZeroMat(cryptoParams, m_ct, s)

	cipherVec := make(crypto.CipherVector, m_ct)
	tmp := make([]float64, m)
	for i := 0; i < s; i++ {
		for j := 0; j < n; j++ {
			ctid := int(j / slots)
			slotid := j % slots

			ct := A[i][ctid]
			ct = crypto.Mask(cryptoParams, ct, slotid, false)
			ct = crypto.InnerSumAll(cryptoParams, crypto.CipherVector{ct})
			for k := range cipherVec {
				cipherVec[k] = ct
			}

			mat.Row(tmp, j, B)
			plainVec, _ := crypto.EncodeFloatVector(cryptoParams, tmp)
			prod := crypto.CPMult(cryptoParams, cipherVec, plainVec)
			out[i] = crypto.CAdd(cryptoParams, out[i], prod)
		}
	}

	return out
}

// Duplicate individual elements in A
func CPMatMult2F(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, B *mat.Dense) crypto.CipherMatrix {
	n, m := B.Dims()
	s := len(A)
	slots := cryptoParams.GetSlots()
	m_ct := int((m-1)/slots) + 1

	if A[0][0].Level() > 2 {
		fmt.Println("Dropping level. Input:", A[0][0].Level())
		A = crypto.DropLevel(cryptoParams, A, 2)
	}
	fmt.Println("A level:", A[0][0].Level())

	accCache := make([]CipherVectorAccV1, s)
	for i := range accCache {
		accCache[i] = NewCipherVectorAccV1(cryptoParams, m_ct)
	}

	cipherVec := make(crypto.CipherVector, m_ct)
	tmp := make([]float64, m)
	var outScale float64
	for j := 0; j < n; j++ {
		mat.Row(tmp, j, B)
		plainVec, _ := crypto.EncodeFloatVector(cryptoParams, tmp)
		ToMontgomeryForm(cryptoParams, plainVec)

		for i := 0; i < s; i++ {
			ctid := int(j / slots)
			slotid := j % slots

			ct := A[i][ctid]
			ct = crypto.Mask(cryptoParams, ct, slotid, false)
			ct = crypto.InnerSumAll(cryptoParams, crypto.CipherVector{ct})
			for k := range cipherVec {
				cipherVec[k] = ct
			}

			outScale = ct.Scale() * plainVec[0].Scale()

			CPMultAccWithoutMRedV1(cryptoParams, cipherVec, plainVec, accCache[i])
		}
	}

	out := make(crypto.CipherMatrix, s)
	for i := range out {
		out[i] = ModularReduceV1(cryptoParams, accCache[i], outScale)
	}
	fmt.Println("Out Scale (log2)", math.Log2(out[0][0].Scale()), "Out Level", out[0][0].Level())

	return out
}

type BlockI8 struct {
	data [][]int8
	r    int
	c    int
}

func NewBlockI8(r, c int) BlockI8 {
	return BlockI8{
		data: make([][]int8, r),
		r:    r,
		c:    c,
	}
}

func (b BlockI8) At(i, j int) float64 {
	return float64(b.data[i][j])
}

func (b BlockI8) Dims() (int, int) {
	return b.r, b.c
}

type BlockF64 mat.Dense

func (b BlockF64) At(i, j int) float64 {
	return b.At(i, j)
}

type Block interface {
	At(int, int) float64
	Dims() (int, int)
}

type BlockVector []Block
type BlockMatrix []BlockVector

func ToBlockMatrix(A *mat.Dense, d int) BlockMatrix {
	r, c := A.Dims()
	br, bc := int((r-1)/d)+1, int((c-1)/d)+1

	out := make(BlockMatrix, br)
	for bi := range out {
		out[bi] = make(BlockVector, bc)
		for bj := range out[bi] {
			i1, i2 := bi*d, Min((bi+1)*d, r)
			j1, j2 := bj*d, Min((bj+1)*d, c)
			out[bi][bj] = A.Slice(i1, i2, j1, j2)
		}
	}

	return out
}

// Return if a diagonal vector exists without extracting elements
func GetDiagBool(X Block, dim int, index int) bool {
	r, c := X.Dims()
	index = Mod(index, dim) // range [0, dim-1]
	return (dim+1-r) <= index || index <= c-1
}

// index 0 is the main diagonal
// max size of Block is dim by dim and index ranges from 0 to dim-1 (mod dim)
// If given diagonal does not overlap with X (matrix might be smaller), returns false
func GetDiag(dst []float64, X Block, dim int, index int) bool {
	r, c := X.Dims()

	index = Mod(index, dim) // range [0, dim-1]

	if (dim+1-r) <= index || index <= c-1 {

		if dst == nil {
			dst = make([]float64, dim)
		} else if len(dst) < c {
			panic("destination array is not large enough")
		}

		i := Mod(-index, dim)
		for j := 0; j < len(dst); j++ {
			if i < r && j < c {
				dst[j] = X.At(i, j)
			} else {
				dst[j] = 0
			}

			i = Mod(i+1, dim)
		}

		return true
	}

	return false
}

func convertToComplex128WithRot(v []float64, nrot int) []complex128 {
	res := make([]complex128, len(v))
	for i, el := range v {
		res[Mod(i+nrot, len(res))] = complex(el, 0)
	}
	return res
}

// Return if a diagonal vector exists without extracting/encoding the vectors
func EncodeDiagBool(X BlockVector, index int, slots int) bool {
	for i := range X {
		if GetDiagBool(X[i], slots, index) {
			return true
		}
	}
	return false
}

// index specifies which diagonal to extract
// applies right-rotation by nrot positions before encoding
func EncodeDiag(cryptoParams *crypto.CryptoParams, X BlockVector, index int, nrot int, level int) (crypto.PlainVector, bool) {
	slots := cryptoParams.GetSlots()

	buf := make([]float64, slots)
	out := make(crypto.PlainVector, len(X))
	anyFlag := false

	for i := range X {
		success := GetDiag(buf, X[i], slots, index)
		if success {
			anyFlag = true
			plaintext := ckks.NewPlaintext(cryptoParams.Params, level, cryptoParams.Params.Scale())
			cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
				encoder.EncodeNTT(plaintext, convertToComplex128WithRot(buf, nrot), cryptoParams.Params.LogSlots())
				return nil
			})
			out[i] = plaintext
		} else {
			out[i] = nil
		}
	}

	return out, anyFlag
}

func EncodeDiagWithEncoder(cryptoParams *crypto.CryptoParams, X BlockVector, index int, nrot int, level int, enc ckks.Encoder) (crypto.PlainVector, bool) {
	slots := cryptoParams.GetSlots()

	buf := make([]float64, slots)
	out := make(crypto.PlainVector, len(X))
	anyFlag := false

	for i := range X {
		success := GetDiag(buf, X[i], slots, index)
		if success {
			anyFlag = true
			plaintext := ckks.NewPlaintext(cryptoParams.Params, level, cryptoParams.Params.Scale())
			enc.EncodeNTT(plaintext, convertToComplex128WithRot(buf, nrot), cryptoParams.Params.LogSlots())
			out[i] = plaintext
		} else {
			out[i] = nil
		}
	}

	return out, anyFlag
}

// Pre-rotate A, all shifts, mult with diagonals from B
func CPMatMult3(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, B *mat.Dense) crypto.CipherMatrix {
	s := len(A)
	slots := cryptoParams.GetSlots()
	blockB := ToBlockMatrix(B, slots)
	fmt.Println("blockB dims:", len(blockB), len(blockB[0]))

	m_ct := len(blockB[0])

	if A[0][0].Level() > 2 {
		fmt.Println("Dropping level. Input:", A[0][0].Level())
		A = crypto.DropLevel(cryptoParams, A, 2)
	}
	fmt.Println("A level:", A[0][0].Level())

	accCache := make([]CipherVectorAccV1, s)
	for i := range accCache {
		accCache[i] = NewCipherVectorAccV1(cryptoParams, m_ct)
	}

	cipherVec := make(crypto.CipherVector, m_ct)

	for bi := range blockB {
		Acol := make(crypto.CipherVector, s)
		for i := range Acol {
			Acol[i] = A[i][bi]
		}

		breakFlag := false
		for shift := 0; shift < slots; shift++ {
			plainVec, flag := EncodeDiag(cryptoParams, blockB[bi], -shift, 0, cryptoParams.Params.MaxLevel())
			if !flag {
				breakFlag = true
				break
			}

			ToMontgomeryForm(cryptoParams, plainVec)

			for i := range Acol {
				for j := range cipherVec {
					cipherVec[j] = Acol[i]
				}
				CPMultAccWithoutMRedV1(cryptoParams, cipherVec, plainVec, accCache[i])
			}

			// Shift left by one
			for i := range Acol {
				Acol[i] = crypto.RotateRight(cryptoParams, Acol[i], -1)
			}
		}

		if breakFlag {
			// Reset Acol
			for i := range Acol {
				Acol[i] = crypto.RotateRight(cryptoParams, A[i][bi], 1)
			}

			for shift := -1; shift > -slots; shift-- {
				plainVec, flag := EncodeDiag(cryptoParams, blockB[bi], -shift, 0, cryptoParams.Params.MaxLevel())
				if !flag {
					break
				}

				ToMontgomeryForm(cryptoParams, plainVec)

				for i := range Acol {
					for j := range cipherVec {
						cipherVec[j] = Acol[i]
					}
					CPMultAccWithoutMRedV1(cryptoParams, cipherVec, plainVec, accCache[i])
				}

				// Shift right by one
				for i := range Acol {
					Acol[i] = crypto.RotateRight(cryptoParams, Acol[i], 1)
				}
			}
		}
	}

	out := make(crypto.CipherMatrix, s)
	outScale := A[0][0].Scale() * cryptoParams.Params.Scale()
	for i := range out {
		out[i] = ModularReduceV1(cryptoParams, accCache[i], outScale)
	}
	return out
}

// Pre-rotate A, baby step giant step version, mult with diagonals from B
func CPMatMult4(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, B *mat.Dense, maxLevel int) crypto.CipherMatrix {

	s := len(A)
	outScale := A[0][0].Scale() * cryptoParams.Params.Scale()
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))
	blockB := ToBlockMatrix(B, slots)
	fmt.Println("slots", slots, "d", d)
	fmt.Println("blockB dims:", len(blockB), len(blockB[0]))

	m_ct := len(blockB[0])

	if A[0][0].Level() > maxLevel {
		fmt.Println("Dropping level. Input:", A[0][0].Level())
		A = crypto.DropLevel(cryptoParams, A, maxLevel)
	}
	fmt.Println("A level:", A[0][0].Level())

	accCache := make([][]CipherVectorAccV2, s)
	for i := range accCache {
		accCache[i] = make([]CipherVectorAccV2, d) // Cache each of the sqrt(slots) groups
		//Initialize later on-the-fly
	}

	for bi := range blockB {

		rotCache := make(crypto.CipherMatrix, s)
		for i := range rotCache {
			rotCache[i] = make(crypto.CipherVector, d)
		}

		for shift := 0; shift < slots; shift++ {
			baby, giant := shift%d, int(shift/d)
			plainVec, flag := EncodeDiag(cryptoParams, blockB[bi], -shift, d*giant, cryptoParams.Params.MaxLevel())
			if !flag {
				continue
			}

			fmt.Println("EncodeDiag", "baby", baby, "giant", giant)

			ToMontgomeryForm(cryptoParams, plainVec)

			for i := range A {

				if rotCache[i][baby] == nil {
					rotCache[i][baby] = crypto.RotateRight(cryptoParams, A[i][bi], -baby)
				}

				cipherVec := make(crypto.CipherVector, m_ct)
				for j := range cipherVec {
					cipherVec[j] = rotCache[i][baby]
				}

				if accCache[i][giant].val == nil {
					accCache[i][giant] = NewCipherVectorAccV2(cryptoParams, m_ct, maxLevel)
				}

				CPMultAccWithoutMRedV2(cipherVec, plainVec, accCache[i][giant])
			}
		}
	}

	out := make(crypto.CipherMatrix, s)
	for i := range out {
		for l := range accCache[i] {
			if accCache[i][l].val != nil {
				cv := ModularReduceV2(cryptoParams, accCache[i][l], outScale)
				if l > 0 { // Giant step alignment
					for j := range cv {
						cv[j] = crypto.RotateRight(cryptoParams, cv[j], -l*d)
					}
				}

				if out[i] == nil {
					out[i] = cv
				} else {
					out[i] = crypto.CAdd(cryptoParams, out[i], cv)
				}
			}
		}
	}
	return out

}

// Cache structure for MatMult4
// First index corresponds to row index of blocks of size slots-by-slots
// Second index of indexMap corresponds to index of diagonals (0..slots-1),
// The table maps a diag index to the encoded PlainVector
// If a given index has empty data, stored element is nil
type PlainMatrixDiagCache [][]crypto.PlainVector

func MatMult4StreamPreprocess(cryptoParams *crypto.CryptoParams, gfs *GenoFileStream, maxLevel int, cacheFilePrefix string) {
	gfs.Reset() // Reset to beginning of file just in case

	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))

	m_ct := ((gfs.NumCols() - 1) / uint64(slots)) + 1
	numBlockRows := ((gfs.NumRows() - 1) / uint64(slots)) + 1
	nproc := runtime.GOMAXPROCS(0)

	log.LLvl1(time.Now().Format(time.RFC3339), "MatMult4StreamPreprocess:", "input", gfs.NumRows(), gfs.NumCols(), "numBlockRows", numBlockRows, "numBlockCols", m_ct)

	for bi := 0; bi < int(numBlockRows); bi++ {

		dcs, flag := NewDiagCacheStream(cryptoParams, cacheFilePrefix, bi, true)
		if flag {
			continue
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "gathering submatrix")

		BSlice := make([]BlockI8, m_ct)
		nr := Min((bi+1)*slots, int(gfs.NumRows())) - bi*slots
		for ri := 0; ri < nr; ri++ {

			// Read one row from file
			row := gfs.NextRow()

			// Add slice to each block matrix
			for bj := range BSlice {
				j1 := bj * slots
				j2 := Min((bj+1)*slots, int(gfs.NumCols()))
				nc := j2 - j1
				if ri == 0 {
					BSlice[bj] = NewBlockI8(nr, nc)
				}
				BSlice[bj].data[ri] = row[j1:j2]
			}
		}

		blockVec := make(BlockVector, m_ct)
		for bj := range blockVec {
			blockVec[bj] = Block(BSlice[bj])
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "finding active diagonals")

		// Pre-collect active baby/giant indices
		babyTable := make([]bool, d)
		giantTable := make([]bool, d)
		shiftTable := make([]bool, slots)
		for shift := 0; shift < slots; shift++ {
			if EncodeDiagBool(blockVec, -shift, slots) {
				baby, giant := shift%d, shift/d
				babyTable[baby] = true
				giantTable[giant] = true
				shiftTable[shift] = true
			}
		}

		dcs.SetIndexTables(babyTable, giantTable)

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "extracting and caching diagonals")

		type dataItem struct {
			plainVec crypto.PlainVector
			shift    int
		}

		jobChannels := make([]chan int, nproc)
		for i := range jobChannels {
			jobChannels[i] = make(chan int, 32)
		}

		diagChannel := make(chan dataItem, 16)

		// Job feeder
		go func() {
			for shift, flag := range shiftTable {
				if flag {
					jobChannels[shift%nproc] <- shift
				}
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		// Data writer
		var writer sync.WaitGroup
		writer.Add(1)
		go func() {
			defer writer.Done()
			for item := range diagChannel {
				dcs.WriteDiag(item.plainVec, uint32(item.shift))
			}
		}()

		// Data encoders
		var encoderGroup sync.WaitGroup
		for thread := 0; thread < nproc; thread++ {
			encoderGroup.Add(1)
			go func(thread int) {
				defer encoderGroup.Done()

				enc := ckks.NewEncoderBig(cryptoParams.Params, cryptoParams.GetPrec())

				for shift := range jobChannels[thread] {
					_, giant := shift%d, shift/d

					plainVec, _ := EncodeDiagWithEncoder(cryptoParams, blockVec, -shift, d*giant, maxLevel, enc)

					ToMontgomeryForm(cryptoParams, plainVec)

					diagChannel <- dataItem{plainVec, shift}
				}
			}(thread)
		}

		encoderGroup.Wait()
		close(diagChannel)

		writer.Wait()
		dcs.Close()
	}

	return
}

func MatMult4StreamCompute(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, maxLevel int, cacheFilePrefix string) crypto.CipherMatrix {
	s := len(A)
	outScale := A[0][0].Scale() * cryptoParams.Params.Scale()
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))
	nproc := runtime.GOMAXPROCS(0)

	numBlockRows := len(A[0])

	log.LLvl1(time.Now().Format(time.RFC3339), "MatMult4StreamCompute")
	if A[0][0].Level() > maxLevel {
		log.LLvl1(time.Now().Format(time.RFC3339), "Dropping level. Input:", A[0][0].Level())
		A = crypto.DropLevel(cryptoParams, A, maxLevel)
	}

	accCache := make([][]CipherVectorAccV2, s)
	accCacheMux := make([][]sync.Mutex, s)
	for i := range accCache {
		accCache[i] = make([]CipherVectorAccV2, d) // Cache each of the sqrt(slots) groups
		accCacheMux[i] = make([]sync.Mutex, d)
	}

	rotCache := make(crypto.CipherMatrix, s)
	for i := range rotCache {
		rotCache[i] = make(crypto.CipherVector, d)
	}

	var m_ct int

	for bi := 0; bi < numBlockRows; bi++ {

		dcs, _ := NewDiagCacheStream(cryptoParams, cacheFilePrefix, bi, false)

		babyTable, giantTable := dcs.GetIndexTables()

		m_ct = int(dcs.vectorLen)

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "generating rotation cache")

		// Dispatcher
		jobChannels := make([]chan int, nproc)
		for i := range jobChannels {
			jobChannels[i] = make(chan int, 32)
		}

		go func() {
			for baby, flag := range babyTable {
				if flag {
					jobChannels[baby%nproc] <- baby
				}
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		// Workers
		var workerGroup sync.WaitGroup
		Aslice := make(crypto.CipherVector, len(A))
		for i := range A {
			Aslice[i] = A[i][bi]
		}
		for thread := 0; thread < nproc; thread++ {
			workerGroup.Add(1)
			go func(thread int) {
				defer workerGroup.Done()

				eva := ckks.NewEvaluator(cryptoParams.Params, ckks.EvaluationKey{Rlk: cryptoParams.Rlk, Rtks: cryptoParams.RotKs})

				for baby := range jobChannels[thread] {
					for i := range A {
						rotCache[i][baby] = crypto.RotateRightWithEvaluator(cryptoParams, Aslice[i], -baby, eva)
					}
				}
			}(thread)
		}
		workerGroup.Wait()

		for giant, flag := range giantTable {
			if flag {
				for i := range A {
					if accCache[i][giant].val == nil {
						accCache[i][giant] = NewCipherVectorAccV2(cryptoParams, int(dcs.vectorLen), maxLevel)
					}
				}
			}
		}

		var wg sync.WaitGroup

		type dataItem struct {
			plainVec crypto.PlainVector
			shift    int
		}

		diagChannels := make([]chan dataItem, nproc)
		for i := range diagChannels {
			diagChannels[i] = make(chan dataItem, 8)
		}

		// Data feeder
		go func() {
			for plainVec, shift := dcs.ReadDiag(); plainVec != nil; plainVec, shift = dcs.ReadDiag() {
				diagChannels[shift%nproc] <- dataItem{plainVec, shift}
			}
			for _, c := range diagChannels {
				close(c)
			}
		}()

		// Data processors
		for thread := 0; thread < nproc; thread++ {
			wg.Add(1)
			go func(thread int) {
				defer wg.Done()
				for item := range diagChannels[thread] {
					plainVec, shift := item.plainVec, item.shift
					baby, giant := shift%d, shift/d
					for i := range A {
						accCacheMux[i][giant].Lock()
						CPMultAccWithoutMRedV2(crypto.CipherVector{rotCache[i][baby]}, plainVec, accCache[i][giant])
						accCacheMux[i][giant].Unlock()
					}
				}
			}(thread)
		}
		wg.Wait()
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Postprocessing accumulators")

	out := crypto.CZeroMat(cryptoParams, m_ct, s)
	for i := range out {
		aggChannel := make(chan crypto.CipherVector, 16)

		jobChannels := make([]chan int, nproc)
		for j := range jobChannels {
			jobChannels[j] = make(chan int, 32)
		}

		go func() {
			for l := range accCache[i] {
				if accCache[i][l].val != nil {
					jobChannels[l%nproc] <- l
				}
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		var wg sync.WaitGroup
		for thread := 0; thread < nproc; thread++ {
			wg.Add(1)
			go func(thread int) {
				defer wg.Done()

				eva := ckks.NewEvaluator(cryptoParams.Params, ckks.EvaluationKey{Rlk: cryptoParams.Rlk, Rtks: cryptoParams.RotKs})

				for l := range jobChannels[thread] {
					cv := ModularReduceV2(cryptoParams, accCache[i][l], outScale)

					if l > 0 { // Giant step alignment
						for j := range cv {
							cv[j] = crypto.RotateRightWithEvaluator(cryptoParams, cv[j], -l*d, eva)
						}
					}

					aggChannel <- cv
				}
			}(thread)
		}

		var aggGroup sync.WaitGroup
		aggGroup.Add(1)
		go func() {
			defer aggGroup.Done()

			eva := ckks.NewEvaluator(cryptoParams.Params, ckks.EvaluationKey{Rlk: cryptoParams.Rlk, Rtks: cryptoParams.RotKs})

			for cv := range aggChannel {
				for j := range cv {
					eva.Add(out[i][j], cv[j], out[i][j])
				}
			}
		}()

		wg.Wait()
		close(aggChannel)
		aggGroup.Wait()
	}

	return out
}

func MatMult4Stream(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, gfs *GenoFileStream, maxLevel int, computeSquaredSum bool, nproc int) (crypto.CipherMatrix, []float64, []float64) {
	gfs.Reset() // Reset to beginning of file just in case

	nrow, ncol := gfs.NumRowsToKeep(), gfs.NumColsToKeep()
	if nproc <= 0 { // If nproc is non-positive, use all cores
		nproc = runtime.GOMAXPROCS(0)
	}

	s := len(A)
	outScale := A[0][0].Scale() * cryptoParams.Params.Scale()
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))

	//blockB := ToBlockMatrix(B, slots)
	//fmt.Println("blockB dims:", len(blockB), len(blockB[0]))
	m_ct := ((ncol - 1) / uint64(slots)) + 1
	numBlockRows := ((nrow - 1) / uint64(slots)) + 1

	if A[0][0].Level() > maxLevel {
		fmt.Println("Dropping level. Input:", A[0][0].Level())
		A = crypto.DropLevel(cryptoParams, A, maxLevel)
	}
	fmt.Println("A level:", A[0][0].Level())

	accCache := make([][]CipherVectorAccV2, s)
	accCacheMux := make([][]sync.Mutex, s)
	for i := range accCache {
		accCache[i] = make([]CipherVectorAccV2, d) // Cache each of the sqrt(slots) groups, initialize later on-the-fly
		accCacheMux[i] = make([]sync.Mutex, d)
	}

	rotCache := make(crypto.CipherMatrix, s)
	for i := range rotCache {
		rotCache[i] = make(crypto.CipherVector, d)
	}

	var sqSum, sum []float64
	if computeSquaredSum {
		sqSum = make([]float64, ncol)
		sum = make([]float64, ncol)
	}

	for bi := 0; bi < int(numBlockRows); bi++ {

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "gathering submatrix")

		BSlice := make([]BlockI8, m_ct)
		nr := Min((bi+1)*slots, int(nrow)) - bi*slots
		for ri := 0; ri < nr; ri++ {

			// Read one row from file
			row := gfs.NextRow()

			// Replace missing with zeros
			for rj := range row {
				if row[rj] < 0 {
					row[rj] = 0
				}

				if computeSquaredSum {
					sqSum[rj] += float64(row[rj] * row[rj])
					sum[rj] += float64(row[rj])
				}
			}

			// Add slice to each block matrix
			for bj := range BSlice {
				j1 := bj * slots
				j2 := Min((bj+1)*slots, int(ncol))
				nc := j2 - j1
				if ri == 0 {
					BSlice[bj] = NewBlockI8(nr, nc)
				}
				BSlice[bj].data[ri] = row[j1:j2]
			}
		}

		blockVec := make(BlockVector, m_ct)
		for bj := range blockVec {
			blockVec[bj] = Block(BSlice[bj])
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "finding active diagonals")

		// Pre-collect active baby/giant indices
		babyTable := make([]bool, d)
		giantTable := make([]bool, d)
		shiftTable := make([]bool, slots)
		for shift := 0; shift < slots; shift++ {
			if EncodeDiagBool(blockVec, -shift, slots) {
				baby, giant := shift%d, shift/d
				babyTable[baby] = true
				giantTable[giant] = true
				shiftTable[shift] = true
			}
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "generating rotation cache")

		log.LLvl1(time.Now().Format(time.RFC3339), "Num procs", nproc)

		// Dispatcher
		jobChannels := make([]chan int, nproc)
		for i := range jobChannels {
			jobChannels[i] = make(chan int, 64)
		}
		go func() {
			index := 0
			for baby, flag := range babyTable {
				if flag {
					jobChannels[index%nproc] <- baby
					index++
				}
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		// Workers
		var workerGroup sync.WaitGroup
		Aslice := make(crypto.CipherVector, len(A))
		for i := range A {
			Aslice[i] = A[i][bi]
		}
		for thread := 0; thread < nproc; thread++ {
			workerGroup.Add(1)
			go func(thread int) {
				defer workerGroup.Done()

				eva := ckks.NewEvaluator(cryptoParams.Params, ckks.EvaluationKey{Rlk: cryptoParams.Rlk, Rtks: cryptoParams.RotKs})

				for baby := range jobChannels[thread] {
					for i := range A {
						rotCache[i][baby] = crypto.RotateRightWithEvaluator(cryptoParams, Aslice[i], -baby, eva)
					}
				}
			}(thread)
		}
		workerGroup.Wait()

		for giant, flag := range giantTable {
			if flag {
				for i := range A {
					if accCache[i][giant].val == nil {
						accCache[i][giant] = NewCipherVectorAccV2(cryptoParams, int(m_ct), maxLevel)
					}
				}
			}
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "extracting and multiplying diagonals")

		// Extract and encode diagonal vectors
		shiftChannels := make([]chan int, nproc)
		for i := range shiftChannels {
			shiftChannels[i] = make(chan int, 128)
		}

		go func() {
			index := 0
			for shift, flag := range shiftTable {
				if flag {
					if (index+1)%1000 == 0 {
						log.LLvl1(index + 1)
					}
					shiftChannels[index%nproc] <- shift
					index++
				}
			}
			for _, c := range shiftChannels {
				close(c)
			}
		}()

		for thread := 0; thread < nproc; thread++ {
			workerGroup.Add(1)
			go func(thread int) {
				defer workerGroup.Done()

				enc := ckks.NewEncoderBig(cryptoParams.Params, cryptoParams.GetPrec())

				for shift := range shiftChannels[thread] {
					baby, giant := shift%d, shift/d

					plainVec, _ := EncodeDiagWithEncoder(cryptoParams, blockVec, -shift, d*giant, maxLevel, enc)

					ToMontgomeryForm(cryptoParams, plainVec)

					for i := range A {
						accCacheMux[i][giant].Lock()
						CPMultAccWithoutMRedV2(crypto.CipherVector{rotCache[i][baby]}, plainVec, accCache[i][giant])
						accCacheMux[i][giant].Unlock()
					}
				}
			}(thread)
		}
		workerGroup.Wait()
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Postprocessing accumulators")

	out := crypto.CZeroMat(cryptoParams, int(m_ct), s)
	for i := range out {
		jobChannels := make([]chan int, nproc)
		for j := range jobChannels {
			jobChannels[j] = make(chan int, 32)
		}

		go func() {
			for l := range accCache[i] {
				if accCache[i][l].val != nil {
					jobChannels[l%nproc] <- l
				}
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		aggChannel := make(chan crypto.CipherVector, 8)

		var wg sync.WaitGroup
		for thread := 0; thread < nproc; thread++ {
			wg.Add(1)
			go func(thread int) {
				defer wg.Done()

				eva := ckks.NewEvaluator(cryptoParams.Params, ckks.EvaluationKey{Rlk: cryptoParams.Rlk, Rtks: cryptoParams.RotKs})

				for l := range jobChannels[thread] {
					cv := ModularReduceV2(cryptoParams, accCache[i][l], outScale)

					if l > 0 { // Giant step alignment
						for j := range cv {
							cv[j] = crypto.RotateRightWithEvaluator(cryptoParams, cv[j], -l*d, eva)
						}
					}

					aggChannel <- cv
				}
			}(thread)
		}

		var aggGroup sync.WaitGroup
		aggGroup.Add(1)
		go func() {
			defer aggGroup.Done()

			eva := ckks.NewEvaluator(cryptoParams.Params, ckks.EvaluationKey{Rlk: cryptoParams.Rlk, Rtks: cryptoParams.RotKs})

			for cv := range aggChannel {
				for j := range cv {
					eva.Add(out[i][j], cv[j], out[i][j])
				}
			}
		}()

		wg.Wait()
		close(aggChannel)
		aggGroup.Wait()
	}

	return out, sum, sqSum
}

func MatMult4TransformB(cryptoParams *crypto.CryptoParams, B *mat.Dense) PlainMatrixDiagCache {
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))
	blockB := ToBlockMatrix(B, slots)

	cache := make(PlainMatrixDiagCache, len(blockB))

	for bi := range blockB {

		cache[bi] = make([]crypto.PlainVector, slots)

		for shift := 0; shift < slots; shift++ {
			giant := int(shift / d)
			plainVec, flag := EncodeDiag(cryptoParams, blockB[bi], -shift, d*giant, cryptoParams.Params.MaxLevel())
			if !flag {
				cache[bi][shift] = nil
			} else {
				ToMontgomeryForm(cryptoParams, plainVec)
				cache[bi][shift] = plainVec
			}
		}
	}

	return cache
}

// Caches a banded matrix B with specified upper and lower bandwidth.
// Assumes that B is at most a (slots x slots) matrix, for `slots` the number of slots in a ciphertext.
func MatMult4TransformBandedB(cryptoParams *crypto.CryptoParams, B *mat.Dense, upperBandwidth, lowerBandwidth int) PlainMatrixDiagCache {
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))
	blockB := ToBlockMatrix(B, slots)

	cache := make(PlainMatrixDiagCache, len(blockB))

	// The band is the diagonals with shift indices in [0, upperBandwidth + 1) and [slots - lowerBandwidth, slots)
	topDiagBound := upperBandwidth + 1 // The +1 is because we also want to include the main diagonal, which is index 0.
	leftDiagBound := slots - lowerBandwidth
	// log.LLvl1("topDiagBound:", topDiagBound, "; leftDiagBound", leftDiagBound)

	for bi := range blockB {

		cache[bi] = make([]crypto.PlainVector, slots)

		for shift := 0; shift < slots; shift++ {
			if topDiagBound <= shift && shift < leftDiagBound {
				// The diagonal this shift corresponds falls outside of the band.
				continue
			}

			giant := int(shift / d)
			plainVec, flag := EncodeDiag(cryptoParams, blockB[bi], -shift, d*giant, cryptoParams.Params.MaxLevel())
			if !flag {
				cache[bi][shift] = nil
			} else {
				ToMontgomeryForm(cryptoParams, plainVec)
				cache[bi][shift] = plainVec
			}
		}
	}

	return cache
}

func CPMatMult4CachedB(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, CachedB PlainMatrixDiagCache) crypto.CipherMatrix {
	s := len(A)
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))
	fmt.Println("slots", slots, "d", d)
	fmt.Println("brows", len(CachedB))

	if A[0][0].Level() > 2 {
		fmt.Println("Dropping level. Input:", A[0][0].Level())
		A = crypto.DropLevel(cryptoParams, A, 2)
	}
	fmt.Println("A level:", A[0][0].Level())

	out := make(crypto.CipherMatrix, s)
	outScale := A[0][0].Scale() * cryptoParams.Params.Scale()

	//rotCount := 1

	for i := range A {

		accCache := make([]CipherVectorAccV1, d) // Cache each of the sqrt(slots) groups

		for bi := range CachedB {

			rotCache := make(crypto.CipherVector, d)

			for shift := 0; shift < slots; shift++ {
				if CachedB[bi][shift] == nil {
					continue
				}

				//fmt.Println("i", i+1, "/", len(A), "bi", bi+1, "/", len(CachedB), "shift", shift)

				baby, giant := shift%d, int(shift/d)
				plainVec := CachedB[bi][shift]

				if rotCache[baby] == nil {
					//fmt.Println("Rot #", rotCount)
					//rotCount = rotCount + 1
					rotCache[baby] = crypto.RotateRight(cryptoParams, A[i][bi], -baby)
				}

				cipherVec := make(crypto.CipherVector, len(plainVec))
				for j := range cipherVec {
					cipherVec[j] = rotCache[baby]
				}

				if accCache[giant] == nil {
					accCache[giant] = NewCipherVectorAccV1(cryptoParams, len(plainVec))
				}

				CPMultAccWithoutMRedV1(cryptoParams, cipherVec, plainVec, accCache[giant])
			}
		}

		for l := range accCache {
			if accCache[l] != nil {

				//fmt.Println("ModReduce accCache: i", i+1, "giant", l)

				cv := ModularReduceV1(cryptoParams, accCache[l], outScale)
				if l > 0 { // Giant step alignment
					for j := range cv {
						//fmt.Println("Rot #", rotCount)
						//rotCount = rotCount + 1
						cv[j] = crypto.RotateRight(cryptoParams, cv[j], -l*d)
					}
				}

				if out[i] == nil {
					out[i] = cv
				} else {
					out[i] = crypto.CAdd(cryptoParams, out[i], cv)
				}
			}
		}
	}
	return out
}

// Generalized to levels >= 2
func CPMatMult4V2CachedB(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, maxLevel int, CachedB PlainMatrixDiagCache) crypto.CipherMatrix {
	s := len(A)
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))
	//fmt.Println("slots", slots, "d", d)
	//fmt.Println("brows", len(CachedB))

	if A[0][0].Level() > maxLevel {
		fmt.Println("Dropping level. Input:", A[0][0].Level())
		A = crypto.DropLevel(cryptoParams, A, maxLevel)
	}
	fmt.Println("CPMatMult4V2CachedB, A level", A[0][0].Level())

	out := make(crypto.CipherMatrix, s)
	outScale := A[0][0].Scale() * cryptoParams.Params.Scale()

	for i := range A {

		accCache := make([]CipherVectorAccV2, d) // Cache each of the sqrt(slots) groups

		for bi := range CachedB {

			rotCache := make(crypto.CipherVector, d)

			for shift := 0; shift < slots; shift++ {
				if CachedB[bi][shift] == nil {
					continue
				}

				baby, giant := shift%d, int(shift/d)
				plainVec := CachedB[bi][shift]

				if rotCache[baby] == nil {
					rotCache[baby] = crypto.RotateRight(cryptoParams, A[i][bi], -baby)
				}

				cipherVec := make(crypto.CipherVector, len(plainVec))
				for j := range cipherVec {
					cipherVec[j] = rotCache[baby]
				}

				if accCache[giant].val == nil {
					accCache[giant] = NewCipherVectorAccV2(cryptoParams, len(plainVec), maxLevel)
				}

				CPMultAccWithoutMRedV2(cipherVec, plainVec, accCache[giant])
			}
		}

		for l := range accCache {
			if accCache[l].val != nil {

				cv := ModularReduceV2(cryptoParams, accCache[l], outScale)
				if l > 0 { // Giant step alignment
					for j := range cv {
						cv[j] = crypto.RotateRight(cryptoParams, cv[j], -l*d)
					}
				}

				if out[i] == nil {
					out[i] = cv
				} else {
					out[i] = crypto.CAdd(cryptoParams, out[i], cv)
				}
			}
		}
	}
	return out
}
