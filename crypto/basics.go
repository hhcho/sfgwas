package crypto

import (
	"math"
	"math/cmplx"
	"time"

	"github.com/ldsec/lattigo/v2/ckks"
	"go.dedis.ch/onet/v3/log"
	"gonum.org/v1/gonum/mat"
)

func EncodeDense(cryptoParams *CryptoParams, vals *mat.Dense) PlainMatrix {
	//col encoded
	r, c := vals.Dims()
	valsEnc := make(PlainMatrix, c)
	tmp := make([]float64, r)
	for i := 0; i < c; i++ {
		valsC := mat.Col(tmp, i, vals)
		valsEnc[i], _ = EncodeFloatVector(cryptoParams, valsC)
	}
	return valsEnc

}

func EncryptDense(cryptoParams *CryptoParams, vals *mat.Dense) CipherMatrix {
	//col encoded
	r, c := vals.Dims()
	valsEnc := make(CipherMatrix, c)
	tmp := make([]float64, r)
	for i := 0; i < c; i++ {
		valsC := mat.Col(tmp, i, vals)
		valsEnc[i], _ = EncryptFloatVector(cryptoParams, valsC)
	}
	return valsEnc
}

//Takes PlaintextMatrix and returns a denseMatrix: if col is true (ie pt is a column encoded matrix, then we transpose)
func PlaintextToDense(cryptoParams *CryptoParams, pt PlainMatrix, ptVecSize int) *mat.Dense {
	vals := make([]float64, len(pt)*ptVecSize)
	for i := range pt {
		tmp := DecodeFloatVector(cryptoParams, pt[i])
		for j := 0; j < ptVecSize; j++ {
			vals[i*ptVecSize+j] = tmp[j]
		}
	}

	denseMat := mat.NewDense(len(pt), ptVecSize, vals)
	return mat.DenseCopyOf(denseMat.T())
}

func EncryptPlaintextMatrix(cryptoParams *CryptoParams, pm PlainMatrix) CipherMatrix {
	cm := make(CipherMatrix, len(pm))
	for c := range pm {
		cm[c] = make(CipherVector, len(pm[c]))
		for i := range pm[c] {
			cryptoParams.WithEncryptor(func(encryptor ckks.Encryptor) error {
				cm[c][i] = encryptor.EncryptNew(pm[c][i])
				return nil
			})
		}
	}
	return cm
}

func DummyBootstrapMatrix(cryptoParams *CryptoParams, Ap CipherMatrix) CipherMatrix {
	for i := range Ap {
		Ap[i].DummyBootstrapping("dummy", cryptoParams)
	}
	return Ap
}

func GlobalToPartyIndex(cryptoParams *CryptoParams, Arowdims []int, col, nparty int) (int, int, int) {
	pid := 0
	ctid := 0
	slotid := col
	for i := 0; i < nparty; i++ {
		if slotid < Arowdims[i] {
			pid = i
			ctid = int(slotid / (cryptoParams.GetSlots()))
			slotid = slotid % (cryptoParams.GetSlots())
			break
		} else {
			slotid = slotid - Arowdims[i]
		}
	}
	return pid, ctid, slotid
}

//DCopyEncrypted returns a shallow? copy of A?
func DCopyEncrypted(A []CipherMatrix) []CipherMatrix {
	Acopy := make([]CipherMatrix, len(A))
	for p := range A {
		Acopy[p] = CopyEncryptedMatrix(A[p])

	}
	return Acopy
}

/* Ciphertext operations */
func GetCTsize(vecsize, ctind, slots int) int {
	if (ctind+1)*slots <= vecsize {
		return slots
	}
	return vecsize % slots
}

// Retain only the first N slots
// TODO is there a way to do this without consuming a level?
func MaskTrunc(cryptoParams *CryptoParams, ct *ckks.Ciphertext, N int) *ckks.Ciphertext {
	if N == cryptoParams.GetSlots() {
		return ct
	}

	m := make([]float64, cryptoParams.GetSlots())
	for i := 0; i < N; i++ {
		m[i] = 1.0
	}

	mask, _ := EncodeFloatVector(cryptoParams, m)
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		ct = eval.MulRelinNew(mask[0], ct)
		error := eval.Rescale(ct, cryptoParams.Params.Scale(), ct)
		return error
	})
	return ct
}

func MaskWithScaling(cryptoParams *CryptoParams, ct *ckks.Ciphertext, ind int, keep bool, scalingFactor float64) *ckks.Ciphertext {
	m := make([]float64, cryptoParams.GetSlots())
	if keep {
		for i := range m {
			if i != ind {
				m[i] = scalingFactor
			}
		}
	} else {
		m[ind] = scalingFactor
	}

	mask, _ := EncodeFloatVector(cryptoParams, m)
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		ct = eval.MulRelinNew(mask[0], ct)
		error := eval.Rescale(ct, cryptoParams.Params.Scale(), ct)
		return error
	})
	return ct
}

func Mask(cryptoParams *CryptoParams, ct *ckks.Ciphertext, index int, keepRest bool) *ckks.Ciphertext {
	if ct == nil {
		return nil
	}
	m := make([]float64, cryptoParams.GetSlots())
	if keepRest {
		for i := range m {
			if i != index {
				m[i] = 1.0
			}
		}
	} else {
		m[index] = 1.0
	}

	mask, _ := EncodeFloatVector(cryptoParams, m)
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		ct = eval.MulRelinNew(mask[0], ct)
		error := eval.Rescale(ct, cryptoParams.Params.Scale(), ct)
		return error
	})
	return ct
}

func Add(cryptoParams *CryptoParams, ct1 *ckks.Ciphertext, ct2 *ckks.Ciphertext) *ckks.Ciphertext {
	var newct *ckks.Ciphertext
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		newct = eval.AddNew(ct1, ct2)
		return nil
	})
	return newct
}

func AddPlain(cryptoParams *CryptoParams, ct1 *ckks.Ciphertext, ct2 *ckks.Plaintext) *ckks.Ciphertext {
	var newct *ckks.Ciphertext
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		newct = eval.AddNew(ct1, ct2)
		return nil
	})
	return newct
}

func AddConst(cryptoParams *CryptoParams, ct *ckks.Ciphertext, constant interface{}) *ckks.Ciphertext {
	var resct *ckks.Ciphertext
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		resct = eval.AddConstNew(ct, constant)
		return nil
	})
	return resct
}

func RotateRightWithEvaluator(cryptoParams *CryptoParams, ct *ckks.Ciphertext, nrot int, eva ckks.Evaluator) *ckks.Ciphertext {
	nrot = Mod(nrot, cryptoParams.GetSlots())
	var newct *ckks.Ciphertext
	if nrot != 0 {
		newct = eva.RotateNew(ct, cryptoParams.GetSlots()-nrot)
	} else {
		newct = ct.CopyNew().Ciphertext()
	}
	return newct
}

func RotateRight(cryptoParams *CryptoParams, ct *ckks.Ciphertext, nrot int) *ckks.Ciphertext {
	nrot = Mod(nrot, cryptoParams.GetSlots())
	var newct *ckks.Ciphertext
	if nrot != 0 {
		cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			newct = eval.RotateNew(ct, cryptoParams.GetSlots()-nrot)
			return nil
		})
	} else {
		newct = ct.CopyNew().Ciphertext()
	}
	return newct
}

func Mult(cryptoParams *CryptoParams, ct1 *ckks.Ciphertext, ct2 *ckks.Ciphertext) *ckks.Ciphertext {
	var res *ckks.Ciphertext
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		res = eval.MulRelinNew(ct1, ct2)
		return nil
	})
	return res
}

//RotateAndAdd computes the inner sum of a Ciphertext
func RotateAndAdd(cryptoParams *CryptoParams, ct *ckks.Ciphertext, size int) *ckks.Ciphertext {
	ctOut := ct.CopyNew().Ciphertext()
	for rotate := 1; rotate < size; rotate *= 2 {
		cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			rt := eval.RotateNew(ctOut, rotate)
			eval.Add(rt, ctOut, ctOut)
			return nil
		})
	}
	return ctOut
}

func Rebalance(cryptoParams *CryptoParams, ct *ckks.Ciphertext) *ckks.Ciphertext {
	if ct == nil {
		return nil
	}
	ctOut := InnerSumAll(cryptoParams, CipherVector{ct})
	ctOut = CMultConst(cryptoParams, CipherVector{ctOut}, 1/float64(cryptoParams.GetSlots()), false)[0]
	return ctOut
}

func InnerSum(cryptoParams *CryptoParams, X CipherVector, Xsize int) *ckks.Ciphertext {
	// Sum all the ciphertexts in vector (i.e vector sum)
	vecsum := X[0].CopyNew().Ciphertext()
	if len(X) == 1 {
		return RotateAndAdd(cryptoParams, X[0], Xsize)
	}

	for i := 1; i < len(X); i++ {
		cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			eval.Add(X[i], vecsum, vecsum)
			return nil
		})
	}

	return RotateAndAdd(cryptoParams, X[0], cryptoParams.GetSlots())
}

func InnerProd(cryptoParams *CryptoParams, X, Y CipherVector) *ckks.Ciphertext {
	return InnerSumAll(cryptoParams, CMult(cryptoParams, X, Y))
}

func InnerSumAll(cryptoParams *CryptoParams, X CipherVector) *ckks.Ciphertext {
	slots := cryptoParams.GetSlots()
	vecsum := X[0].CopyNew().Ciphertext()

	if len(X) > 1 {
		for i := 1; i < len(X); i++ {
			cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
				eval.Add(X[i], vecsum, vecsum)
				return nil
			})
		}
	}

	return RotateAndAdd(cryptoParams, vecsum, slots)
}

// RotateAndPlace rotates ct to the right by shift = place
func RotateAndPlace(cryptoParams *CryptoParams, ct *ckks.Ciphertext, size, place int, duplicate bool) *ckks.Ciphertext {
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		if !duplicate {
			if !(place == 0) {
				//rotate right
				ct = eval.RotateNew(ct, cryptoParams.GetSlots()-place)
			}

		} else {
			maxnrot := int(math.Ceil(math.Log2(float64(size))))
			rotate := 1
			for j := 0; j < maxnrot; j++ {
				tmp := eval.RotateNew(ct, cryptoParams.GetSlots()-rotate)
				ct = eval.AddNew(ct, tmp)
				rotate = rotate * 2
			}
		}
		return nil
	})
	return ct

}

// AggregateSumMask sums across parties
func AggregateSumMask(cryptoParams *CryptoParams, vals []*ckks.Ciphertext) *ckks.Ciphertext {
	aggsum := vals[0].CopyNew().Ciphertext()
	if len(vals) == 1 {
		return aggsum
	}

	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := 1; i < len(vals); i++ {
			eval.Add(vals[i], aggsum, aggsum)
		}
		return nil
	})
	aggsum = Mask(cryptoParams, aggsum, 0, false)
	return aggsum
}

func AggregateVec(cryptoParams *CryptoParams, vecs []CipherVector) CipherVector {
	res := CopyEncryptedVector(vecs[0])
	for i := 1; i < len(vecs); i++ {
		res = CAdd(cryptoParams, vecs[i], res)
	}
	return res
}

func AggregateMat(cryptoParams *CryptoParams, mats []CipherMatrix) CipherMatrix {
	res := CopyEncryptedMatrix(mats[0])

	for i := 1; i < len(mats); i++ {
		for j := range mats[0] {
			res[j] = CAdd(cryptoParams, mats[i][j], res[j])
		}
	}
	return res
}

// SqSum evaluates the sum(X^2) across the vector
func SqSum(cryptoParams *CryptoParams, X CipherVector) *ckks.Ciphertext {
	X2 := CMult(cryptoParams, X, X)
	return InnerSumAll(cryptoParams, X2)
}

/* Ciphervector operations */
func Zero(cryptoParams *CryptoParams) *ckks.Ciphertext {
	a := []float64{0.0}
	tmp, _ := EncryptFloatVector(cryptoParams, a)
	return tmp[0]
}

func CZeros(cryptoParams *CryptoParams, n int) CipherVector {
	cv := make(CipherVector, n)
	a := []float64{0.0}
	for i := range cv {
		tmp, _ := EncryptFloatVector(cryptoParams, a)
		cv[i] = tmp[0]
	}
	return cv
}

//column based
func CZeroMat(cryptoParams *CryptoParams, nrows, ncols int) CipherMatrix {
	cm := make(CipherMatrix, ncols)
	for r := range cm {
		cm[r] = CZeros(cryptoParams, nrows)
	}
	return cm
}

func CMult(cryptoParams *CryptoParams, X CipherVector, Y CipherVector) CipherVector {
	xlen := len(X)
	ylen := len(Y)
	res := make(CipherVector, Max(len(X), len(Y)))
	if xlen == 1 {
		cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			for i := range Y {
				res[i] = eval.MulRelinNew(Y[i], X[0])
				error := eval.Rescale(res[i], cryptoParams.Params.Scale(), res[i])
				if error != nil {
					return error
				}
			}
			return nil
		})
		return res
	}
	if ylen == 1 {
		cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			for i := range X {
				res[i] = eval.MulRelinNew(X[i], Y[0])
				error := eval.Rescale(res[i], cryptoParams.Params.Scale(), res[i])
				if error != nil {
					return error
				}
			}
			return nil
		})
		return res
	}
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range X {
			res[i] = eval.MulRelinNew(X[i], Y[i])
			error := eval.Rescale(res[i], cryptoParams.Params.Scale(), res[i])
			if error != nil {
				return error
			}
		}
		return nil
	})
	return res
}

func CPMult(cryptoParams *CryptoParams, X CipherVector, Y PlainVector) CipherVector {
	xlen := len(X)
	ylen := len(Y)
	res := make(CipherVector, Max(len(X), len(Y)))
	if xlen == 1 {
		cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			for i := range Y {
				res[i] = eval.MulRelinNew(Y[i], X[0])
				error := eval.Rescale(res[i], cryptoParams.Params.Scale(), res[i])
				if error != nil {
					return error
				}
			}
			return nil
		})
		return res
	}
	if ylen == 1 {
		cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			for i := range X {
				res[i] = eval.MulRelinNew(X[i], Y[0])
				error := eval.Rescale(res[i], cryptoParams.Params.Scale(), res[i])
				if error != nil {
					return error
				}
			}
			return nil
		})
		return res
	}
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range X {
			res[i] = eval.MulRelinNew(X[i], Y[i])
			error := eval.Rescale(res[i], cryptoParams.Params.Scale(), res[i])
			if error != nil {
				return error
			}
		}
		return nil
	})
	return res
}

func CMultConstMat(cryptoParams *CryptoParams, X CipherMatrix, constant interface{}, inPlace bool) (res CipherMatrix) {
	res = make(CipherMatrix, len(X))
	for i := range res {
		res[i] = CMultConst(cryptoParams, X[i], constant, inPlace)
	}
	return res
}

func CMultConst(cryptoParams *CryptoParams, X CipherVector, constant interface{}, inPlace bool) (res CipherVector) {
	if inPlace {
		res = X
	} else {
		res = CopyEncryptedVector(X)
	}
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range X {
			if inPlace {
				eval.MultByConst(res[i], constant, res[i])
			} else {
				res[i] = eval.MultByConstNew(res[i], constant)
			}
		}
		return nil
	})
	return res
}

func FlattenLevels(cryptoParams *CryptoParams, X CipherMatrix) (CipherMatrix, int) {
	minLevel := math.MaxInt32
	notFlat := false
	for i := range X {
		for j := range X[i] {
			if X[i][j].Level() != minLevel {
				minLevel = Min(X[i][j].Level(), minLevel)
				notFlat = true
			}
		}
	}

	if !notFlat {
		return X, minLevel
	}

	return DropLevel(cryptoParams, X, minLevel), minLevel
}

func CMultConstRescale(cryptoParams *CryptoParams, X CipherVector, constant interface{}, inPlace bool) CipherVector {
	var res CipherVector
	if inPlace {
		res = X
	} else {
		res = CopyEncryptedVector(X)
	}
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range X {
			res[i] = eval.MultByConstNew(res[i], constant)
			error := eval.Rescale(res[i], cryptoParams.Params.Scale(), res[i])
			if error != nil {
				return error
			}
		}
		return nil
	})
	return res
}

func CMultScalar(cryptoParams *CryptoParams, X CipherVector, ct *ckks.Ciphertext) CipherVector {
	res := CopyEncryptedVector(X)
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range X {
			res[i] = eval.MulRelinNew(res[i], ct)
			error := eval.Rescale(res[i], cryptoParams.Params.Scale(), res[i])
			if error != nil {
				return error
			}
		}
		return nil
	})
	return res
}

func CAdd(cryptoParams *CryptoParams, X CipherVector, Y CipherVector) CipherVector {
	res := make(CipherVector, len(X)) //equal num of ciphertexts
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := 0; i < Max(len(Y), len(X)); i++ {
			//check level
			res[i] = eval.AddNew(X[i], Y[i])
		}
		return nil
	})
	return res
}

func CSub(cryptoParams *CryptoParams, X CipherVector, Y CipherVector) CipherVector {
	res := make(CipherVector, len(X))
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range Y {
			res[i] = eval.SubNew(X[i], Y[i])
		}
		return nil
	})
	return res

}

func CPAdd(cryptoParams *CryptoParams, X CipherVector, Y PlainVector) CipherVector {
	res := make(CipherVector, len(X))
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range Y {
			res[i] = eval.AddNew(X[i], Y[i])
		}
		return nil
	})
	return res

}

func CAddConst(cryptoParams *CryptoParams, X CipherVector, constant interface{}) CipherVector {
	res := CopyEncryptedVector(X)
	for i := range X {
		res[i] = AddConst(cryptoParams, X[i], constant)
	}
	return res

}

func ChebyApproximation(cryptoParams *CryptoParams, X CipherVector, cheby *ckks.ChebyshevInterpolation) CipherVector {
	res := CopyEncryptedVector(X)
	err := cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		var err error
		for i := range X {
			res[i], err = eval.EvaluateCheby(res[i], cheby, res[i].Scale())
			return err
		}
		return nil
	})

	if err != nil {
		log.Fatal(err)
	}
	return res
}

func CInverse(cryptoParams *CryptoParams, X CipherVector, intv IntervalApprox) CipherVector {
	y := make(CipherVector, len(X))
	err := cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range X {
			y[i] = eval.InverseNew(X[i], intv.Iter)
		}
		return nil
	})
	if err != nil {
		log.Fatal(err)
	}
	return y
}

func CRescale(cryptoParams *CryptoParams, X CipherVector) CipherVector {
	err := cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := range X {
			eval.Rescale(X[i], cryptoParams.Params.Scale(), X[i])
		}
		return nil
	})
	if err != nil {
		log.Fatal(err)
	}
	return X
}

func LevelTest(ciphers CipherVector, cryptoParams *CryptoParams, needed int, serverID, name string) CipherVector {
	if ciphers[0].Level() <= needed {
		log.LLvl1(time.Now().Format(time.RFC3339), "Dummy Bootstrapping called")
		ciphers.DummyBootstrapping(serverID, cryptoParams)

	}
	return ciphers
}

func LevelTestMatrix(ciphers CipherMatrix, cryptoParams *CryptoParams, needed int, serverID, name string) CipherMatrix {
	if ciphers[0][0].Level() <= needed {
		log.LLvl1(time.Now().Format(time.RFC3339), "Dummy Bootstrapping required for function ", name)
		for i := range ciphers {
			ciphers[i].DummyBootstrapping(serverID, cryptoParams)
		}
	}
	return ciphers
}

/*TODO: LEGACY CODE NEED TO REMOVE*/

// InvSqrtApprox computes an encrypted approximated version of the inverse function
func InvSqrtApprox(cryptoParams *CryptoParams, ctIn *ckks.Ciphertext, intv IntervalApprox) *ckks.Ciphertext {
	cheby := ckks.Approximate(invSqrt, complex(intv.A, 0), complex(intv.B, 0), intv.Degree)

	var res *ckks.Ciphertext
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		var err error
		res, err = eval.EvaluateCheby(ctIn, cheby, cryptoParams.Params.Scale())
		return err
	})
	return res
}

// SqrtApprox computes an encrypted approximated version of the inverse function
func SqrtApprox(cryptoParams *CryptoParams, ctIn *ckks.Ciphertext, intv IntervalApprox) *ckks.Ciphertext {
	cheby := ckks.Approximate(Sqrt, complex(intv.A, 0), complex(intv.B, 0), intv.Degree)

	var res *ckks.Ciphertext
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		var err error
		res, err = eval.EvaluateCheby(ctIn, cheby, cryptoParams.Params.Scale())
		return err
	})
	return res
}

func invSqrt(x complex128) complex128 {
	return 1 / cmplx.Sqrt(x)
}

func Sqrt(x complex128) complex128 {
	return cmplx.Sqrt(x)
}

// InvApprox computes an encrypted approximated version of the inverse function
func InvApprox(cryptoParams *CryptoParams, ctIn *ckks.Ciphertext, intv IntervalApprox) *ckks.Ciphertext {
	cheby := ckks.Approximate(inv, complex(intv.A, 0), complex(intv.B, 0), intv.Degree)

	var approx *ckks.Ciphertext
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		var err error
		approx, err = eval.EvaluateCheby(ctIn, cheby, cryptoParams.Params.Scale())
		return err
	})
	return approx
}

func inv(x complex128) complex128 {
	return 1 / x
}

func ConcatCipherMatrix(mats []CipherMatrix) CipherMatrix {
	var nrow int
	for i := range mats {
		if mats[i] != nil {
			nrow = len(mats[i])
			break
		}
	}

	ncol := 0
	for i := range mats {
		if mats[i] != nil {
			ncol += len(mats[i][0])
		}
	}

	out := make(CipherMatrix, nrow)
	for i := range out {
		out[i] = make(CipherVector, ncol)
		shift := 0
		for j := range mats {
			if mats[j] != nil {
				for c := range mats[j][i] {
					out[i][shift+c] = mats[j][i][c]
				}
				shift += len(mats[j][i])
			}
		}
	}

	return out
}

func DropLevel(cryptoParams *CryptoParams, A CipherMatrix, outLevel int) CipherMatrix {
	out := make(CipherMatrix, len(A))
	for i := range out {
		out[i] = make(CipherVector, len(A[i]))
		for j := range out[i] {
			cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
				if A[i][j].Level() > outLevel {
					out[i][j] = eval.DropLevelNew(A[i][j], A[i][j].Level()-outLevel)
				} else if A[i][j].Level() == outLevel {
					out[i][j] = A[i][j].CopyNew().Ciphertext()
				} else {
					log.Fatalf("DropLevel: requested level", outLevel, "when input is", A[i][j].Level())
				}
				return nil
			})
		}
	}
	return out
}

// Returns the complex conjugate of a ciphervector.
func ComplexConjugate(cryptoParams *CryptoParams, X CipherVector) CipherVector {
	res := make(CipherVector, len(X)) // equal num of ciphertexts
	cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
		for i := 0; i < len(X); i++ {
			// check level
			res[i] = eval.ConjugateNew(X[i])
		}
		return nil
	})
	return res
}

// Returns the real component of a CipherVector.
func CReal(cps *CryptoParams, Z CipherVector) (real CipherVector) {
	conjZ := ComplexConjugate(cps, Z)
	twiceReal := CAdd(cps, Z, conjZ)
	CMultConst(cps, twiceReal, 0.5, true)
	real = CRescale(cps, twiceReal)
	return
}
