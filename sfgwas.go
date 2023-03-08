package main

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strconv"

	"github.com/BurntSushi/toml"
	"github.com/hhcho/sfgwas-private/gwas"
	"github.com/hhcho/sfgwas-private/crypto"
	"github.com/raulk/go-watchdog"

	"github.com/ldsec/lattigo/v2/ckks"
)

// Expects a party ID provided as an environment variable;
// e.g., run "PID=1 go run sfgwas.go"
var PID, PID_ERR = strconv.Atoi(os.Getenv("PID"))

// Default config path
var CONFIG_PATH = "config/"

func main() {
	// RunGWAS()
	playground()
	
}

func InitProtocol(configPath string) *gwas.ProtocolInfo {
	config := new(gwas.Config)

	// Import global parameters
	if _, err := toml.DecodeFile(filepath.Join(configPath, "configGlobal.toml"), config); err != nil {
		fmt.Println(err)
		return nil
	}

	// Import local parameters
	if _, err := toml.DecodeFile(filepath.Join(configPath, fmt.Sprintf("configLocal.Party%d.toml", PID)), config); err != nil {
		fmt.Println(err)
		return nil
	}

	// Create cache/output directories
	if err := os.MkdirAll(config.CacheDir, 0755); err != nil {
		panic(err)
	}
	if err := os.MkdirAll(config.OutDir, 0755); err != nil {
		panic(err)
	}

	// Set max number of threads
	runtime.GOMAXPROCS(config.LocalNumThreads)

	return gwas.InitializeGWASProtocol(config, PID, false)
}

func RunGWAS() {
	if PID_ERR != nil {
		panic(PID_ERR)
	}

	// Initialize protocol
	prot := InitProtocol(CONFIG_PATH)

	// Invoke memory manager
	err, stopFn := watchdog.HeapDriven(prot.GetConfig().MemoryLimit, 40, watchdog.NewAdaptivePolicy(0.5))
	if err != nil {
		panic(err)
	}
	defer stopFn()

	// Run protocol
	prot.GWAS()

	prot.SyncAndTerminate(true)
}

// function for just playing around with the different available tools, manipulating things in ciphertext space, etc.
func playground() {
	precision := uint(8) // 8 and 32 seem to work well, while 16 for some reason seems to accumulate fairly large error (10^-3 ish)
	params := crypto.NewCryptoParamsForNetwork(ckks.DefaultParams[ckks.PN13QP218], 1, precision)
	selfParams := params[0] // the params for network suppose multiple parties, but we just need one party, so we just take the first set of params

	// get rotation keys
	selfParams.SetRotKeys(crypto.GenerateRotKeys(selfParams.GetSlots(), 20, true))

	// size := 4
	// plaintextVector1 := make([]float64, size)
	// plaintextVector2 := make([]float64, size)
	// for i := 0; i < size; i++ {
	// 	plaintextVector1[i] = float64(i*i)
	// 	plaintextVector2[i] = float64(-i)
	// }
	// fmt.Println(plaintextVector1)
	// fmt.Println(plaintextVector2)
	// encryptedVector1, _ := crypto.EncryptFloatVector(selfParams, plaintextVector1)
	// encryptedVector2, _ := crypto.EncryptFloatVector(selfParams, plaintextVector2)
	// encSum := crypto.InnerSumAll(selfParams, encryptedVector1)
	// decSum := crypto.DecryptFloat(selfParams, encSum)
	// fmt.Println(decSum)

	// encryptedResult := crypto.CAdd(selfParams, encryptedVector1, encryptedVector2) // CAdd standing for "Ciphertext Addition"
	// DecryptFloatVector both "Decrypts" (puts back into polynomial encoding) and "Decodes" (returns to message space)
	// note that, perhaps slightly confusingly, we call the polynomial encoding space the "plaintexts"
	// plaintextResult := crypto.DecryptFloatVector(selfParams, encryptedResult, size)
	// fmt.Println(plaintextResult)

	// multiply the matrices
	// 1 2 5      -2  0
	// 3 4 6       4 -1
	//             6  5
	// expected result:
	// 36 23
	// 46 26
	// A := [][]float64{
	// 	{1, 2, 5},
	// 	{3, 4, 6},
	// }
	// B := [][]float64{
	// 	{-2, 0},
	// 	{4, -1},
	// 	{6, 5},
	// }
	A := [][]float64{
		{1, 2},
		{2, 4},
	}
	B := [][]float64{
		{-4, 0},
		{2, -1},
	}
	fmt.Println(B)
	B = transposeMatrix(B)
	fmt.Println(B)

	// flattenedA := flattenMatrix(A)
	// flattenedB := flattenMatrix(B)

	// fmt.Println(flattenedA)
	// fmt.Println(flattenedB)

	// encA, _ := crypto.EncryptFloatVector(selfParams, flattenedA)
	// fmt.Printf("encrypted %d elements\n", i)
	// fmt.Println(len(encA))
	// encB, _ := crypto.EncryptFloatVector(selfParams, flattenedB)

	encA := encryptMatrixByRows(selfParams, A)
	encB := encryptMatrixByRows(selfParams, B)

	encResult := CMatrixMultiply2(selfParams, encA, encB, 2)
	plainResult := crypto.DecryptFloatVector(selfParams, encResult, 4)
	fmt.Println(plainResult)
	matrixResult := unflattenMatrix(plainResult, 2)
	fmt.Println(matrixResult)

	// encA = crypto.RotateAndAdd(selfParams, encA, 6)

	// doubleA := crypto.CipherVector{encA[0], encA[0]}
	// doubleASquared := crypto.CMult(selfParams, doubleA, doubleA)
	// fmt.Println(crypto.DecryptFloatVector(selfParams, doubleASquared, 12))

	// tmp2 := crypto.MaskTrunc(selfParams, encA[0], 3)

	// tmp := crypto.CMult(selfParams, encA, crypto.CipherVector{tmp2})
	// fmt.Println(crypto.DecryptFloatVector(selfParams, tmp, 6))


	// decA := crypto.DecryptFloatVector(selfParams, crypto.CipherVector{encA[0]}, 6)
	// fmt.Println(decA)

	// fmt.Println(crypto.DecryptFloatVector(selfParams,
	// 	crypto.CipherVector{crypto.RotateAndAdd(selfParams, encA[0], 6)},
	// 	8))

	// crucially the granularity on RotateAndAdd is powers of 2: that is to say, you can only request that it sum the indices 0 through 2^k. Any upper index greater than 2^k returns the same value as 2^k+1. (Up to minor fluctuations induced by the "empty" slots having extremely small values in them.) Check the algorithm for RotateAndAdd for details.
	// the upshot is that if you want to say sum just the first n indices of a vector, then you should mask the vector first (with MaskTrunc)
	
	// decUniqueA := crypto.DecryptFloat(selfParams, encA[0])
	// fmt.Println(decUniqueA)

}

func encryptMatrixByRows(cps *crypto.CryptoParams, A [][]float64) crypto.CipherVector {
	result := make(crypto.CipherVector, 0)
	for i := 0; i < len(A); i++ {
		row, _ := crypto.EncryptFloatVector(cps, A[i])
		result = append(result, row[0])
	}
	return result
}

// computes and returns a Ciphertext with the value of the sum of the first n slots in a given ciphertext
// TODO alternatively you can actually call crypto.RotateAndAdd, repackage its output to be a CipherVector and decrypt it as a Float Vector to look at the slots beyond the first
// by construction of the algorithm, the ith slot in the Ciphertext output of RotateAndAdd is the sum excluding terms 0 through i, so it's possible to recover the partial sums of the first few terms this way
func truncatedRotateAndAdd(cps *crypto.CryptoParams, ct *ckks.Ciphertext, n int) *ckks.Ciphertext {
	return crypto.RotateAndAdd(cps, crypto.MaskTrunc(cps, ct, n), n)
}

// flatten an n x m matrix A and return it as a single slice (of length n*m)
// the matrix is represented in a "row packed" way, that is to say that the rows are kept whole and laid end-to-end (while the columns get spread out)
func flattenMatrix(A [][]float64) []float64 {
	n, m := len(A), len(A[0])
	flattenedA := make([]float64, n*m)
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			flattenedA[(i*m) + j] = A[i][j]
		}
	}
	return flattenedA
}

// unflatten a flattened matrix into a new matrix with `n` rows
func unflattenMatrix(v []float64, n int) [][]float64 {
	// first, make sure it's possible to flatten to this dimension at all
	l := len(v)
	if l % n != 0 {
		errorString := fmt.Sprintf("length %d vector cannot be unflattened into a matrix of %d rows", l, n)
		panic(errorString)
	}
	m := l / n

	matrix := make([][]float64, n)
	for i := 0; i < n; i++ {
		matrix[i] = v[(i*m):(i+1)*m]
	}
	return matrix
}

// transpose a matrix A
func transposeMatrix(A [][]float64) [][]float64 {
	n, m := len(A), len(A[0])
	transpose := make([][]float64, m)
	for j := 0; j < m; j++ {
		row := make([]float64, n)
		for i := 0; i < n; i++ {
			row[i] = A[i][j]
		}
		transpose[j] = row
	}
	return transpose
}

// multiply encrypted matrix A : n x k with plaintext matrix B: k x m
func cpMatrixMultiply(params *crypto.CryptoParams, A crypto.CipherMatrix, B [][]float64) crypto.CipherMatrix {
	n := len(A)
	k, m := len(B), len(B[0])
	
	// first, we collect Plaintext encodings of the columns of B
	plaintextCols := make(crypto.PlainMatrix, m)
	for j := 0; j < m; j++ {
		col := make([]float64, k)
		for i := 0; i < k; i++ {
			col[i] = B[i][j]
		}
		// fmt.Print("col is ", col, "\n")
		plainCol, _ := crypto.EncodeFloatVector(params, col);
		plaintextCols[j] = plainCol
	}
	
	// now use the columns of B to compute the result
	result := make(crypto.CipherMatrix, n)
	for i := 0; i < n; i++ {
		row := make(crypto.CipherVector, m)
		for j := 0; j < m; j++ {
			prod := crypto.CPMult(params, A[i], plaintextCols[j])
			entry := crypto.InnerSum(params, prod, k)
			row[j] = entry.CopyNew().Ciphertext()

			fmt.Printf("\njust did entry %d,%d\n", i, j)
			fmt.Print("which multiplied ", crypto.DecryptFloatVector(params, A[i], k), " into ", crypto.DecodeFloatVector(params, plaintextCols[j])[:k],"\n")
			fmt.Print("entry is ", crypto.DecryptFloat(params, row[j]), "\n")
		}
		fmt.Print("row ", i, " is ", crypto.DecryptFloatVector(params, row, 1), "\n")
		// fmt.Print("entry is ", crypto.DecryptFloat(params, row[1]))
		result[i] = row
	}
	fmt.Print("result is ", crypto.DecryptFloatMatrix(params, result, m), "\n")
	return result
}

func encryptMatrixAsVector(cps *crypto.CryptoParams, A [][]float64) crypto.CipherVector {
	n := len(A)
	output := make(crypto.CipherVector, n)
	for i := 0; i < n; i++ {
		enc, _ := crypto.EncryptFloatVector(cps, A[i])
		if len(enc) > 1 {
			// in other words, that we had to pack multiple Ciphertexts
			panic("Encrypted matrix has rows too wide to store in individual Ciphertexts")
		}
		// otherwise we know that there is just one Ciphertext in `enc`, which therefore contains all the values
		output[i] = enc[0]
	}
	return output
}

// given two n x n matrices A and B, computes their matrix product AB, assuming that A is row-encoded (each row is its own *ckks.Ciphertext) and B is column-encoded (each column is its own *ckks.Ciphertext)
func CMatrixMultiply2(cps *crypto.CryptoParams, A, B crypto.CipherVector, n int) crypto.CipherVector {
	// ensure the matrices really are n x n
	if n != len(A) || n != len(B) {
		errorString := fmt.Sprintf("Expected %dx%d matrices, got %dx%d and %dx%d", n, n, len(A), n, n, len(B))
		panic(errorString)
	}

	AB := make(crypto.CipherVector, n)
	// just to be sure, zero out each of the entries
	for i := 0; i < n; i++ {
		AB[i] = crypto.Zero(cps)
	}
	
	diagonals := make(crypto.CipherMatrix, n) // semantically it may make more sense to write []CipherVector

	for k := 0; k < n; k++ {
		diagonal := generateDiagonal(cps, A, B, n, k)
		diagonals[k] = diagonal
	}
	// mask the diagonals in the pertinent way to construct the correct outputs for AB
	for i := 0; i < n; i++ {
		// we want to populate row AB[i] with the pertinent values, which are dispersed over the diagonal vectors
		for k := 0; k < n; k++ {
			// for row i, the kth diagonal contributes its value at column (i + k) % n
			offset := (i + k) % n
			onehot := makeOneHot(cps, n, offset)
			maskedDiagonal := crypto.Mult(cps, diagonals[i][k], onehot)

			// the maskedDiagonal has the diagonal value in its 0th index, so we need to rotate it over to the correct position before continuing
			maskedDiagonal = crypto.RotateRight(cps, maskedDiagonal, offset)
			AB[i] = crypto.Add(cps, AB[i], maskedDiagonal)
			// TODO finish
		}
	}

	return AB
}

// helper function for CMatrixMultiply:
// computes the diagonal of matrix values in AB that are offset +k from the diagonal. For example, k = 0 computes the actual diagonal, k=1 computes the diagonal of values at coordinates (0, 1) (1, 2), (2, 3) ... wrapping around to include (n-1, 0). (Using 0 indexing on the matrix.)
// in general, the kth diagonal is those values whose coordinates are (i, (i+k)%n)
// the returned CipherVector is such that its ith Ciphertext represents the value the diagonal contains in the ith row
func generateDiagonal(cps *crypto.CryptoParams, A, B crypto.CipherVector, n, k int) crypto.CipherVector {
	if !(0 <= k && k < n) {
		errorString := fmt.Sprintf("k must be in range [0..n-1], got k=%d and n=%d", k, n)
		panic(errorString)
	}
	// TODO if we change to a column convention it might make it way easier to pull the values out of the diagonal CipherVectors later when reconstructing the rows of AB
	// we use the convention that the ith index of the output is the diagonal entry whose _row_ is index i
	result := make(crypto.CipherVector, n)
	for i := 0; i < n; i++ {
		row := A[i]
		col := B[(i + k) % n]
		// take the dot product of this row/column
		product := crypto.Mult(cps, row, col)
		result[i] = crypto.RotateAndAdd(cps, product, n)
	}
	return result
}

// creates a Ciphertext representing the encryption of the length-n one hot vector with a single 1 at index k
func makeOneHot(cps *crypto.CryptoParams, n, k int) *ckks.Ciphertext {
	vector := make([]float64, n)
	vector[k] = 1
	result, _ := crypto.EncryptFloatVector(cps, vector)
	return result[0]
}

// computes and returns the 
// func maskDiagonal(cps *crypto.CryptoParams, A crypto.CipherVector, n, k int) crypto.CipherVector {

// }

// create and return a new CipherVector whose entries are the entries of `A` cycled by `k`. (I deliberately use the word "cycle" to avoid "rotate")
// for example, cycling a matrix with rows r1, r2, r3, r4 by 3 returns the matrix with rows r4, r1, r2, r3
func cycleRows(A crypto.CipherVector, k int) crypto.CipherVector {
	n := len(A)
	cycledA := make(crypto.CipherVector, n)
	for i := 0; i < n; i++ {
		offsetRow := (i + k) % n
		cycledA[i] = A[offsetRow].CopyNew().Ciphertext()
	}
	return cycledA
}

// given CipherVector `A` representing an array of length-m vectors, create and return a new CipherVector whose entries represent, for each vector in `A`, the summation of the m entries of that vector
func innerSumRows(cps *crypto.CryptoParams, A crypto.CipherVector, m int) crypto.CipherVector {
	n := len(A)
	innerSumA := make(crypto.CipherVector, n)
	for i := 0; i < n; i++ {
		innerSumA[i] = crypto.RotateAndAdd(cps, A[i], m)
	}
	return innerSumA
}


// high level algorithm idea:
// from A : n x n and B : n x n, to compute AB: n x n
// 1. row pack A, i.e. the ciphertext
//      c_A  = [r1-- |r2-- |...|rn-- |...]
//    where r_i is the ith row of A
// 2. column pack B, i.e. the ciphertext
//      c_B^T = [c1-- |c2-- |...|cn-- |...]
//    where c_j is the jth column of B
// 3. for reasons that will become clearer later, reduplicate the encoded part of c_B^T, i.e.
//      c_B^T = [c1-- |c2-- |...|cn-- |c1-- |c2-- |...|cn-- |...]
// 4. for i in [0..n-1]:
//    4a. rotate c_B^T leftward by n places (so that it's rotated i*n relative to original)
//    4b. compute d_i as the pointwise product of c_A and (shifted) c_B^T
//    4c. some introspection will reveal that

// multiplies row-encoded (n x n) matrix A with column-encoded (n x n) matrix B, returning row-encoded matrix AB
// func cpMatrixMultiply3(cps *crypto.CryptoParams, A, B *ckks.Ciphertext, n int) *ckks.Ciphertext {
	
// 	crypto.RotateRight()

// 	for i := 0; i < n; i++ {

// 	}
// }

// takes a row-encoded matrix `A` and 
// TODO need to think of a better name
// func sumDotProducts(cps *crypto.CryptoParams, A *ckks.Ciphertext, n, k int) *ckks.Ciphertext {
// 	// my first implementation of this will be the extremely non-optimized version where you just iteratively add each position together, accumulating the values into one slot


// 	// one silver lining of doing it this way is that it lets us put the result just into the correct offset position quite easily
// 	// of course, if we had done it properly via binary merge, then 


// 	// small optimization to be made: if we also pass in the "offset" k, then we could possibly
// 	// also determine which position
// }

// masks the given row-encoded matrix `A` and
// TODO think of better name
// func maskDotProducts(cps *crypto.CryptoParams, A *ckks.Ciphertext, n, k int) *ckks.Ciphertext {
	
// }

