package crypto

import (
	"bytes"
	"encoding"
	"encoding/gob"
	"errors"
	"fmt"
	"math"
	"runtime"
	"sort"
	"strconv"
	"sync"
	"time"

	"github.com/ldsec/lattigo/v2/ring"
	"go.dedis.ch/onet/v3/log"

	"github.com/ldsec/lattigo/v2/ckks"
	libunlynx "github.com/ldsec/unlynx/lib"
	"go.dedis.ch/onet/v3/simul/monitor"
)

type IntervalApprox struct {
	A          float64
	B          float64
	Degree     int
	Iter       int
	InverseNew bool
	ScaleDowns []float64
}

// CipherVector is a slice of Ciphertexts
type CipherVector []*ckks.Ciphertext

// CipherMatrix is a slice of slice of Ciphertexts
type CipherMatrix []CipherVector

// PlainVector is a slice of Plaintexts
type PlainVector []*ckks.Plaintext

// PlainMatrix is a slice of slice of Plaintexts
type PlainMatrix []PlainVector

// CryptoParams aggregates all ckks scheme information
type CryptoParams struct {
	Sk          *ckks.SecretKey
	AggregateSk *ckks.SecretKey
	Pk          *ckks.PublicKey
	Rlk         *ckks.RelinearizationKey
	RotKs       *ckks.RotationKeySet
	Params      *ckks.Parameters

	encoders   chan ckks.Encoder
	encryptors chan ckks.Encryptor
	decryptors chan ckks.Decryptor
	evaluators chan ckks.Evaluator

	numThreads int
	prec       uint
	BootLimit  int
}

// CryptoParamsForNetwork stores all crypto info to save to file
type CryptoParamsForNetwork struct {
	params      *ckks.Parameters
	sk          []*ckks.SecretKey
	aggregateSk *ckks.SecretKey
	pk          *ckks.PublicKey
	rlk         *ckks.EvaluationKey
	rotKs       *ckks.RotationKeySet
}

var SideRight = true
var SideLeft = false

// RotationType defines how much we should rotate and in which direction
type RotationType struct {
	Value int
	Side  bool
}

// CKKSParamsForTests are _unsecure_ and fast parameters
//var CKKSParamsForTests = ckks.NewParametersFromLogModuli(8, 7, 1<<30, ckks.LogModuli{LogQi: []uint64{36, 30, 30, 30, 30, 30, 30, 30, 30, 30}, LogPi: []uint64{32, 32, 32}}, 3.2)

// #------------------------------------#
// #------------ INIT ------------------#
// #------------------------------------#

// NewCryptoParams initializes CryptoParams with the given values
func NewCryptoParams(params *ckks.Parameters, sk, aggregateSk *ckks.SecretKey, pk *ckks.PublicKey, rlk *ckks.RelinearizationKey, prec uint, numThreads int) *CryptoParams {
	evaluators := make(chan ckks.Evaluator, numThreads)
	for i := 0; i < numThreads; i++ {
		evalKey := ckks.EvaluationKey{
			Rlk:  rlk,
			Rtks: nil,
		}
		evaluators <- ckks.NewEvaluator(params, evalKey)
	}

	encoders := make(chan ckks.Encoder, numThreads)
	var wg sync.WaitGroup
	for i := 0; i < numThreads; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			encoders <- ckks.NewEncoderBig(params, prec)
		}()
	}
	wg.Wait()

	encryptors := make(chan ckks.Encryptor, numThreads)
	for i := 0; i < numThreads; i++ {
		encryptors <- ckks.NewEncryptorFromPk(params, pk)
	}

	decryptors := make(chan ckks.Decryptor, numThreads)
	for i := 0; i < numThreads; i++ {
		decryptors <- ckks.NewDecryptor(params, aggregateSk)
	}

	return &CryptoParams{
		Params:      params,
		Sk:          sk,
		AggregateSk: aggregateSk,
		Pk:          pk,
		Rlk:         rlk,

		encoders:   encoders,
		encryptors: encryptors,
		decryptors: decryptors,
		evaluators: evaluators,

		numThreads: numThreads,
		prec:       prec,
	}
}

// SetDecryptors sets the decryptors in the CryptoParams object
func (cp *CryptoParams) SetDecryptors(params *ckks.Parameters, sk *ckks.SecretKey) {
	decryptors := make(chan ckks.Decryptor, cp.numThreads)
	for i := 0; i < cp.numThreads; i++ {
		decryptors <- ckks.NewDecryptor(params, sk)
	}
	cp.decryptors = decryptors
}

func (cp *CryptoParams) SetEvaluators(params *ckks.Parameters, rlk *ckks.RelinearizationKey, rtks *ckks.RotationKeySet) {
	evaluators := make(chan ckks.Evaluator, cp.numThreads)
	for i := 0; i < cp.numThreads; i++ {
		evalKey := ckks.EvaluationKey{
			Rlk:  rlk,
			Rtks: rtks,
		}
		evaluators <- ckks.NewEvaluator(params, evalKey)
	}
	cp.evaluators = evaluators
}

// NewCryptoParamsForNetwork initializes a set of nbrNodes CryptoParams each containing: keys, encoder, encryptor, decryptor, etc.
func NewCryptoParamsForNetwork(params *ckks.Parameters, nbrNodes int, prec uint) []*CryptoParams {
	kgen := ckks.NewKeyGenerator(params)

	aggregateSk := ckks.NewSecretKey(params)
	skList := make([]*ckks.SecretKey, nbrNodes)
	rq, _ := ring.NewRing(params.N(), append(params.Qi(), params.Pi()...))

	for i := 0; i < nbrNodes; i++ {
		skList[i] = kgen.GenSecretKey()
		rq.Add(aggregateSk.Value, skList[i].Value, aggregateSk.Value)
	}
	pk := kgen.GenPublicKey(aggregateSk)

	ret := make([]*CryptoParams, nbrNodes)
	for i := range ret {
		rlk := kgen.GenRelinearizationKey(aggregateSk)
		ret[i] = NewCryptoParams(params, skList[i], aggregateSk, pk, rlk, prec, runtime.GOMAXPROCS(0))
	}

	return ret
}

// SetRotKeys sets/adds new rotation keys
func (cp *CryptoParams) SetRotKeys(nbrRot []RotationType) []int {
	kgen := ckks.NewKeyGenerator(cp.Params)
	ks := make([]int, 0)
	for i := range nbrRot {
		fmt.Println("GenRot", i, "/", len(nbrRot))

		var rotation int
		if nbrRot[i].Side == SideRight {
			rotation = cp.GetSlots() - nbrRot[i].Value
		} else {
			rotation = nbrRot[i].Value
		}

		// check if key is already in list
		in := false
		for j := range ks {
			if ks[j] == rotation {
				in = true
				break
			}
		}
		if !in {
			ks = append(ks, rotation)
		}
	}

	cp.RotKs = kgen.GenRotationKeysForRotations(ks, false, cp.AggregateSk)
	cp.SetEvaluators(cp.Params, cp.Rlk, cp.RotKs)
	sort.Ints(ks)
	return ks
}

//
//// GenRot generates a left or right rotation of a ciphertext
//func GenRot(cryptoParams *CryptoParams, rotation int, side bool, rotKeys *ckks.RotationKeys) *ckks.RotationKeys {
//	kgen := ckks.NewKeyGenerator(cryptoParams.Params)
//	if rotKeys == nil {
//		rotKeys = ckks.NewRotationKeys()
//	}
//
//	if side == SideRight {
//		kgen.GenRot(ckks.RotationLeft, cryptoParams.AggregateSk, uint64(cryptoParams.GetSlots()-rotation), rotKeys)
//	} else {
//		kgen.GenRot(ckks.RotationLeft, cryptoParams.AggregateSk, uint64(rotation), rotKeys)
//	}
//	return rotKeys
//}

// Generate rotKeys for power of two shifts up to # of slots
// and for every shift up to smallDim
func GenerateRotKeys(slots int, smallDim int, babyFlag bool) []RotationType {
	rotations := make([]RotationType, 0)

	l := slots
	l = FindClosestPow2(l)

	rot := 1
	for i := 0; i < int(math.Ceil(math.Log2(float64(l)))); i++ {
		rotations = append(rotations, RotationType{
			Value: rot,
			Side:  false,
		})
		rotations = append(rotations, RotationType{
			Value: rot,
			Side:  true,
		})
		rot = rot * 2
	}

	//for baby-step giant-step rotations
	if babyFlag {
		rootl := int(math.Ceil(math.Sqrt(float64(slots))))
		for i := 1; i < rootl; i++ {
			rotations = append(rotations, RotationType{
				Value: i,
				Side:  false,
			})
			rotations = append(rotations, RotationType{
				Value: i * rootl,
				Side:  false,
			})
		}
	}

	// for moving the innersum value to its new position
	for i := 1; i < smallDim; i++ {
		rotations = append(rotations, RotationType{
			Value: i,
			Side:  true,
		})
	}

	return rotations
}

func (cp *CryptoParams) GetPrec() uint {
	return cp.prec
}

// GetSlots gets the number of encodable slots (N/2)
func (cp *CryptoParams) GetSlots() int {
	return cp.Params.Slots()
}

// WithEncoder run the given function with an encoder
func (cp *CryptoParams) WithEncoder(act func(ckks.Encoder) error) error {
	encoder := <-cp.encoders
	err := act(encoder)
	cp.encoders <- encoder
	return err
}

// WithEncryptor run the given function with an encryptor
func (cp *CryptoParams) WithEncryptor(act func(ckks.Encryptor) error) error {
	encryptor := <-cp.encryptors
	err := act(encryptor)
	cp.encryptors <- encryptor
	return err
}

// WithDecryptor run the given function with a decryptor
func (cp *CryptoParams) WithDecryptor(act func(act ckks.Decryptor) error) error {
	decryptor := <-cp.decryptors
	err := act(decryptor)
	cp.decryptors <- decryptor
	return err
}

// WithEvaluator run the given function with an evaluator
func (cp *CryptoParams) WithEvaluator(act func(ckks.Evaluator) error) error {
	eval := <-cp.evaluators
	err := act(eval)
	cp.evaluators <- eval
	return err
}

// #------------------------------------#
// #------------ ENCRYPTION ------------#
// #------------------------------------#

// EncryptFloat encrypts one float64 value.
func EncryptFloat(cryptoParams *CryptoParams, num float64) *ckks.Ciphertext {
	slots := cryptoParams.GetSlots()
	plaintext := ckks.NewPlaintext(cryptoParams.Params, cryptoParams.Params.MaxLevel(), cryptoParams.Params.Scale())

	cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
		encoder.Encode(plaintext, ConvertVectorFloat64ToComplex(PadVector([]float64{num}, slots)), cryptoParams.Params.LogSlots())
		return nil
	})

	var ciphertext *ckks.Ciphertext
	cryptoParams.WithEncryptor(func(encryptor ckks.Encryptor) error {
		ciphertext = encryptor.EncryptNew(plaintext)
		return nil
	})
	return ciphertext
}

// EncryptFloatVector encrypts a slice of float64 values in multiple batched ciphertexts.
// and return the number of encrypted elements.
func EncryptFloatVector(cryptoParams *CryptoParams, f []float64) (CipherVector, int) {
	nbrMaxCoef := cryptoParams.GetSlots()
	length := len(f)

	cipherArr := make(CipherVector, 0)
	elementsEncrypted := 0
	for elementsEncrypted < length {
		start := elementsEncrypted
		end := elementsEncrypted + nbrMaxCoef

		if end > length {
			end = length
		}
		plaintext := ckks.NewPlaintext(cryptoParams.Params, cryptoParams.Params.MaxLevel(), cryptoParams.Params.Scale())
		// pad to 0s
		cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
			encoder.Encode(plaintext, ConvertVectorFloat64ToComplex(PadVector(f[start:end], nbrMaxCoef)), cryptoParams.Params.LogSlots())
			return nil
		})
		var cipher *ckks.Ciphertext
		cryptoParams.WithEncryptor(func(encryptor ckks.Encryptor) error {
			cipher = encryptor.EncryptNew(plaintext)
			return nil
		})
		cipherArr = append(cipherArr, cipher)
		elementsEncrypted = elementsEncrypted + (end - start)
	}
	return cipherArr, elementsEncrypted
}

func EncryptFloatVectorWithScale(cryptoParams *CryptoParams, f []float64, scale float64) (CipherVector, int) {
	nbrMaxCoef := cryptoParams.GetSlots()
	length := len(f)

	cipherArr := make(CipherVector, 0)
	elementsEncrypted := 0
	for elementsEncrypted < length {
		start := elementsEncrypted
		end := elementsEncrypted + nbrMaxCoef

		if end > length {
			end = length
		}
		plaintext := ckks.NewPlaintext(cryptoParams.Params, cryptoParams.Params.MaxLevel(), scale)
		// pad to 0s
		cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
			encoder.Encode(plaintext, ConvertVectorFloat64ToComplex(PadVector(f[start:end], nbrMaxCoef)), cryptoParams.Params.LogSlots())
			return nil
		})
		var cipher *ckks.Ciphertext
		cryptoParams.WithEncryptor(func(encryptor ckks.Encryptor) error {
			cipher = encryptor.EncryptNew(plaintext)
			return nil
		})
		cipherArr = append(cipherArr, cipher)
		elementsEncrypted = elementsEncrypted + (end - start)
	}
	return cipherArr, elementsEncrypted
}

// EncryptFloatMatrixRow encrypts a matrix of float64 to multiple packed ciphertexts.
// For this specific matrix encryption each row is encrypted in a set of ciphertexts.
func EncryptFloatMatrixRow(cryptoParams *CryptoParams, matrix [][]float64) (CipherMatrix, int, int, error) {
	nbrRows := len(matrix)
	d := len(matrix[0])

	matrixEnc := make([]CipherVector, 0)
	for _, row := range matrix {
		if d != len(row) {
			return nil, 0, 0, errors.New("this is not a matrix (expected " + strconv.FormatInt(int64(d), 10) +
				" dimensions but got " + strconv.FormatInt(int64(len(row)), 10))
		}
		rowEnc, _ := EncryptFloatVector(cryptoParams, row)
		matrixEnc = append(matrixEnc, rowEnc)
	}
	return matrixEnc, nbrRows, d, nil
}

// EncodeFloatVector encodes a slice of float64 values in multiple batched plaintext (ready to be encrypted).
// It also returns the number of encoded elements.
func EncodeFloatVectorWithScale(cryptoParams *CryptoParams, f []float64, scale float64) (PlainVector, int) {
	nbrMaxCoef := cryptoParams.GetSlots()
	length := len(f)

	plainArr := make(PlainVector, 0)
	elementsEncoded := 0
	for elementsEncoded < length {
		start := elementsEncoded
		end := elementsEncoded + nbrMaxCoef

		if end > length {
			end = length
		}
		plaintext := ckks.NewPlaintext(cryptoParams.Params, cryptoParams.Params.MaxLevel(), scale)
		cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
			encoder.EncodeNTT(plaintext, ConvertVectorFloat64ToComplex(PadVector(f[start:end], nbrMaxCoef)), cryptoParams.Params.LogSlots())
			return nil
		})
		plainArr = append(plainArr, plaintext)
		elementsEncoded = elementsEncoded + (end - start)
	}
	return plainArr, elementsEncoded
}

// EncodeFloatVector encodes a slice of float64 values in multiple batched plaintext (ready to be encrypted).
// It also returns the number of encoded elements.
func EncodeFloatVector(cryptoParams *CryptoParams, f []float64) (PlainVector, int) {
	nbrMaxCoef := cryptoParams.GetSlots()
	length := len(f)

	plainArr := make(PlainVector, 0)
	elementsEncoded := 0
	for elementsEncoded < length {
		start := elementsEncoded
		end := elementsEncoded + nbrMaxCoef

		if end > length {
			end = length
		}
		plaintext := ckks.NewPlaintext(cryptoParams.Params, cryptoParams.Params.MaxLevel(), cryptoParams.Params.Scale())
		cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
			encoder.EncodeNTT(plaintext, ConvertVectorFloat64ToComplex(PadVector(f[start:end], nbrMaxCoef)), cryptoParams.Params.LogSlots())
			return nil
		})
		plainArr = append(plainArr, plaintext)
		elementsEncoded = elementsEncoded + (end - start)
	}
	return plainArr, elementsEncoded
}

// EncodeFloatMatrixRow encodes a matrix of float64 to multiple packed plaintexts.
// For this specific matrix encoding each row is encoded in a set of plaintexts.
func EncodeFloatMatrixRow(cryptoParams *CryptoParams, matrix [][]float64) (PlainMatrix, int, int, error) {
	nbrRows := len(matrix)
	d := len(matrix[0])

	matrixEnc := make(PlainMatrix, 0)
	for _, row := range matrix {
		if d != len(row) {
			return nil, 0, 0, errors.New("this is not a matrix (expected " + strconv.FormatInt(int64(d), 10) +
				" dimensions but got " + strconv.FormatInt(int64(len(row)), 10))
		}

		rowEnc, _ := EncodeFloatVector(cryptoParams, row)
		matrixEnc = append(matrixEnc, rowEnc)
	}
	return matrixEnc, nbrRows, d, nil
}

// #------------------------------------#
// #------------ DECRYPTION ------------#
// #------------------------------------#

// DecryptFloat decrypts a ciphertext with one float64 value.
func DecryptFloat(cryptoParams *CryptoParams, cipher *ckks.Ciphertext) float64 {
	var ret float64
	var plaintext *ckks.Plaintext

	cryptoParams.WithDecryptor(func(decryptor ckks.Decryptor) error {
		plaintext = decryptor.DecryptNew(cipher)
		return nil
	})
	cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
		ret = real(encoder.Decode(plaintext, cryptoParams.Params.LogSlots())[0])
		return nil
	})

	return ret
}

// DecryptMultipleFloat decrypts a ciphertext with multiple float64 values.
// If nbrEl<=0 it decrypts everything without caring about the number of encrypted values.
// If nbrEl>0 the function returns N elements from the decryption.
func DecryptMultipleFloat(cryptoParams *CryptoParams, cipher *ckks.Ciphertext, nbrEl int) []float64 {
	var plaintext *ckks.Plaintext

	cryptoParams.WithDecryptor(func(decryptor ckks.Decryptor) error {
		plaintext = decryptor.DecryptNew(cipher)
		return nil
	})

	var val []complex128
	cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
		val = encoder.Decode(plaintext, cryptoParams.Params.LogSlots())
		return nil
	})
	dataDecrypted := ConvertVectorComplexToFloat64(val)
	if nbrEl <= 0 {
		return dataDecrypted
	}
	return dataDecrypted[:nbrEl]
}

// DecryptFloatVector decrypts multiple batched ciphertexts with N float64 values and appends
// all data into one single float vector.
// If nbrEl<=0 it decrypts everything without caring about the number of encrypted values.
// If nbrEl>0 the function returns N elements from the decryption.
func DecryptFloatVector(cryptoParams *CryptoParams, fEnc CipherVector, N int) []float64 {
	var plaintext *ckks.Plaintext

	dataDecrypted := make([]float64, 0)
	for _, cipher := range fEnc {
		cryptoParams.WithDecryptor(func(decryptor ckks.Decryptor) error {
			plaintext = decryptor.DecryptNew(cipher)
			return nil
		})
		var val []complex128
		cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
			val = encoder.Decode(plaintext, cryptoParams.Params.LogSlots())
			return nil
		})
		dataDecrypted = append(dataDecrypted, ConvertVectorComplexToFloat64(val)...)
	}

	if N <= 0 {
		return dataDecrypted
	}
	return dataDecrypted[:N]
}

// DecryptFloatMatrix decrypts a matrix (kind of) of multiple packed ciphertexts.
// For this specific matrix decryption each row is encrypted in a set of ciphertexts.
// d is the number of column values
func DecryptFloatMatrix(cryptoParams *CryptoParams, matrixEnc []CipherVector, d int) [][]float64 {
	matrix := make([][]float64, 0)
	for _, rowEnc := range matrixEnc {
		row := DecryptFloatVector(cryptoParams, rowEnc, d)
		matrix = append(matrix, row)
	}
	return matrix
}

// DecodeFloatVector decodes a slice of plaintext values in multiple float64 values.
func DecodeFloatVector(cryptoParams *CryptoParams, fEncoded PlainVector) []float64 {
	dataDecoded := make([]float64, 0)
	for _, plaintext := range fEncoded {
		var val []complex128
		cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
			val = encoder.Decode(plaintext, cryptoParams.Params.LogSlots())
			return nil
		})
		dataDecoded = append(dataDecoded, ConvertVectorComplexToFloat64(val)...)
	}
	return dataDecoded
}

// DecodeFloatVector2 decodes a slice of plaintext values in multiple float64 values, along with their imaginary components.
func DecodeFloatVector2(cryptoParams *CryptoParams, fEncoded PlainVector) [][]complex128 {
	dataDecoded := make([][]complex128, 0)
	for _, plaintext := range fEncoded {
		var val []complex128
		cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
			val = encoder.Decode(plaintext, cryptoParams.Params.LogSlots())
			return nil
		})
		dataDecoded = append(dataDecoded, val)
	}
	return dataDecoded
}

// #------------------------------------#
// #------------ MARSHALL --------------#
// #------------------------------------#

func (cm *CipherMatrix) MarshalBinary() ([]byte, [][]int, error) {
	b := make([]byte, 0)
	ctSizes := make([][]int, len(*cm))
	for i, v := range *cm {
		tmp, n, err := v.MarshalBinary()
		ctSizes[i] = n
		if err != nil {
			return nil, nil, err
		}
		b = append(b, tmp...)
	}

	return b, ctSizes, nil

}

func (cm *CipherMatrix) UnmarshalBinary(cryptoParams *CryptoParams, f []byte, ctSizes [][]int) error {
	*cm = make([]CipherVector, len(ctSizes))

	start := 0
	for i := range ctSizes {
		rowSize := 0
		for j := range ctSizes[i] {
			rowSize += ctSizes[i][j]
		}
		end := start + rowSize
		cv := make(CipherVector, 0)
		// log.LLvl1(time.Now().Format(time.RFC3339), "vector: ", i)
		err := cv.UnmarshalBinary(cryptoParams, f[start:end], ctSizes[i])
		if err != nil {
			return err
		}
		start = end
		(*cm)[i] = cv
	}
	return nil
}

func (cv *CipherVector) MarshalBinary() ([]byte, []int, error) {
	data := make([]byte, 0)
	ctSizes := make([]int, 0)
	for _, ct := range *cv {
		b, err := ct.MarshalBinary()
		if err != nil {
			return nil, nil, err
		}

		// log.LLvl1(time.Now().Format(time.RFC3339), "ct level byte: ", i, b[:8]) //first 8 bytes

		data = append(data, b...)
		ctSizes = append(ctSizes, len(b))
	}
	return data, ctSizes, nil
}

// UnmarshalBinary -> CipherVector: converts an array of bytes to an array of ciphertexts.
func (cv *CipherVector) UnmarshalBinary(cryptoParams *CryptoParams, f []byte, fSizes []int) error {
	*cv = make(CipherVector, len(fSizes))

	start := 0
	for i := 0; i < len(fSizes); i++ {
		ct := ckks.NewCiphertext(cryptoParams.Params, 1, cryptoParams.Params.MaxLevel(), cryptoParams.Params.Scale())

		// log.LLvl1(time.Now().Format(time.RFC3339), "ct level byte: ", i, f[start:start+8]) //first 8 bytes

		if err := ct.UnmarshalBinary(f[start : start+fSizes[i]]); err != nil {
			return err
		}
		(*cv)[i] = ct
		start += fSizes[i]
	}
	return nil
}

type cryptoParamsMarshalable struct {
	Params      *ckks.Parameters
	Sk          []*ckks.SecretKey
	AggregateSk *ckks.SecretKey
	Pk          *ckks.PublicKey
	Rlk         *ckks.RelinearizationKey
	RotKs       *ckks.RotationKeySet
}

// #------------------------------------#
// #-------------- COPY ----------------#
// #------------------------------------#

// CopyEncryptedVector does a copy of an array of ciphertexts to a newly created array
func CopyEncryptedVector(src CipherVector) CipherVector {
	dest := make(CipherVector, len(src))
	for i := 0; i < len(src); i++ {
		if src[i] == nil {
			log.LLvl1(time.Now().Format(time.RFC3339), "nil pointer", i)
		}
		dest[i] = (*src[i]).CopyNew().Ciphertext()
	}
	return dest
}

// CopyEncryptedMatrix does a copy of a matrix of ciphertexts to a newly created array
func CopyEncryptedMatrix(src []CipherVector) []CipherVector {
	dest := make([]CipherVector, len(src))
	for i := 0; i < len(src); i++ {
		dest[i] = CopyEncryptedVector(src[i])
	}
	return dest
}

// DummyBootstrapping mimics the bootstrapping
func (cv *CipherVector) DummyBootstrapping(serverID string, cryptoParams *CryptoParams) CipherVector {
	var bootstrappingTime *monitor.TimeMeasure
	if serverID != "" {
		bootstrappingTime = libunlynx.StartTimer(serverID + "_DummyBootstrapping")
	}

	for i, ct := range *cv {
		decryptedValues := DecryptMultipleFloat(cryptoParams, ct, 0)
		encryptedValues, _ := EncryptFloatVector(cryptoParams, decryptedValues)
		(*cv)[i] = encryptedValues[0]
	}

	if serverID != "" {
		libunlynx.EndTimer(bootstrappingTime)
	}
	return *cv
}

/*******/
/*EDITS*/
/*******/

var _ encoding.BinaryMarshaler = new(CryptoParams)
var _ encoding.BinaryUnmarshaler = new(CryptoParams)

//MarshalBinary for minimal cryptoParams-keys + params
func (cp *CryptoParams) MarshalBinary() ([]byte, error) {
	var ret bytes.Buffer
	encoder := gob.NewEncoder(&ret)

	if cp.Params == nil {
		log.LLvl1(time.Now().Format(time.RFC3339), "encoding params is nil")

	} else if cp.Sk == nil {
		log.LLvl1(time.Now().Format(time.RFC3339), "encoding Sk is nil")

	} else if cp.AggregateSk == nil {
		log.LLvl1(time.Now().Format(time.RFC3339), "encoding aggregate sk is nil")

	} else if cp.Rlk == nil {
		log.LLvl1(time.Now().Format(time.RFC3339), "encoding Rlk is nil")
	} else if cp.RotKs == nil {
		log.LLvl1(time.Now().Format(time.RFC3339), "encoding Rotks are nil")
	}

	err := encoder.Encode(cryptoParamsMarshalable{
		Params:      cp.Params,
		Sk:          []*ckks.SecretKey{cp.Sk},
		AggregateSk: cp.AggregateSk,
		Pk:          cp.Pk,
		Rlk:         cp.Rlk,
		RotKs:       cp.RotKs,
	})
	if err != nil {
		return nil, fmt.Errorf("encode minimal crypto params: %v", err)
	}

	return ret.Bytes(), nil
}

func (cp *CryptoParams) UnmarshalBinary(data []byte) error {
	decoder := gob.NewDecoder(bytes.NewBuffer(data))

	decodeParams := new(cryptoParamsMarshalable)
	if err := decoder.Decode(decodeParams); err != nil {
		return fmt.Errorf("decode minimal crypto params: %v", err)
	}

	cp.Params = decodeParams.Params
	cp.Sk = decodeParams.Sk[0]
	cp.AggregateSk = decodeParams.AggregateSk
	cp.Pk = decodeParams.Pk
	cp.Rlk = decodeParams.Rlk
	cp.RotKs = decodeParams.RotKs
	return nil
}

// ConvertVectorFloat64ToComplex converts an array of floats to complex
func ConvertVectorFloat64ToComplex(v []float64) []complex128 {
	res := make([]complex128, len(v))
	for i, el := range v {
		res[i] = complex(el, 0)
	}
	return res
}

// ConvertVectorComplexToFloat64 converts an array of complex to float
func ConvertVectorComplexToFloat64(v []complex128) []float64 {
	res := make([]float64, len(v))
	for i, el := range v {
		res[i] = real(el)
	}
	return res
}

// PadVector pads the vector with 0's before encoding/encryption
func PadVector(v []float64, slots int) []float64 {
	toAdd := make([]float64, slots-len(v))
	return append(v, toAdd...)
}

// FindClosestPow2 finds the closest power of 2 bigger than a number n
func FindClosestPow2(n int) int {
	// find closest power of two
	var bigPower2 int
	for bigPower2 = 1; bigPower2 < n; bigPower2 *= 2 {
	}
	return bigPower2
}
