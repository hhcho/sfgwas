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
	I := [][]float64{
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1},
	}
	A := [][]float64{
		{24, 19, 5, 18},
		{22, 34, 9, -23},
		{8, 0, -8, 0},
		{0, -4, -18, 14},
	}
	B := [][]float64{
		{3, -4, 23, 10},
		{2, -1, -8, 0},
		{-5, 0, 4, -2},
		{5, -9, 1, 3},
	}
	B = transposeMatrix(B)
	n := len(B)

	encA := encryptMatrixByRows(selfParams, A)
	encB := encryptMatrixByRows(selfParams, B)
	encI := encryptMatrixByRows(selfParams, I)

	// sanity check
	fmt.Println(decryptMatrixVectorByRows(selfParams, encA, n))
	fmt.Println(decryptMatrixVectorByRows(selfParams, encB, n))

	encResult := CMatrixMultiply(selfParams, encA, encB, n)
	plainResult := decryptMatrixVectorByRows(selfParams, encResult, n)
	fmt.Println(plainResult)
	fmt.Println(decryptMatrixVectorByRows(selfParams, CMatrixMultiply(selfParams, encA, encI, n), n))
	fmt.Println(decryptMatrixVectorByRows(selfParams, CMatrixMultiply(selfParams, encI, encB, n), n))

	// crucially the granularity on RotateAndAdd is powers of 2: that is to say, you can only request that it sum the indices 0 through 2^k. Any upper index greater than 2^k returns the same value as 2^k+1. (Up to minor fluctuations induced by the "empty" slots having extremely small values in them.) Check the algorithm for RotateAndAdd for details.
	// the upshot is that if you want to say sum just the first n indices of a vector, then you should mask the vector first (with MaskTrunc)

}

func encryptMatrixByRows(cps *crypto.CryptoParams, A [][]float64) crypto.CipherVector {
	result := make(crypto.CipherVector, 0)
	for i := 0; i < len(A); i++ {
		row, _ := crypto.EncryptFloatVector(cps, A[i])
		result = append(result, row[0])
	}
	return result
}

// decrypts a CipherVector interpreted as a series of stacked rows, each row represented with a single Ciphertext containing `n` elements
func decryptMatrixVectorByRows(cps *crypto.CryptoParams, A crypto.CipherVector, n int) [][]float64 {
	result := make([][]float64, 0)
	for i := 0; i < len(A); i++ {
		row := crypto.DecryptFloatVector(cps, crypto.CipherVector{A[i]}, n)
		result = append(result, row)
	}
	return result
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

// given two n x n matrices A and B, computes their matrix product AB, assuming that A is row-encoded (each row is its own *ckks.Ciphertext) and B is column-encoded (each column is its own *ckks.Ciphertext)
// TODO generalize to non-square matrices. Diagonal helper function will be a little harder to think about but should work still.
// worst comes to worst, you can always snap to the largest common square value and zero out the other entries
// TODO use CipherMatrices instead of CipherVectors to let the rows be of arbitrary width
func CMatrixMultiply(cps *crypto.CryptoParams, A, B crypto.CipherVector, n int) crypto.CipherVector {
	// ensure the matrices really are n x n
	if n != len(A) || n != len(B) {
		errorString := fmt.Sprintf("Expected %dx%d matrices, got %dx%d and %dx%d", n, n, len(A), n, n, len(B))
		panic(errorString)
	}
	
	diagonals := make(crypto.CipherMatrix, n) // semantically it may make more sense to write []CipherVector
	
	for k := 0; k < n; k++ {
		diagonal := generateDiagonal(cps, A, B, n, k)
		// fmt.Print("diagonal: ")
		// fmt.Println(decryptMatrixVectorByRows(cps, diagonal, 1))
		diagonals[k] = diagonal
	}

	AB := make(crypto.CipherVector, n)
	// mask the diagonals in the pertinent way to construct the correct outputs for AB
	for i := 0; i < n; i++ {
		// we want to populate row AB[i] with the pertinent values, which are dispersed over the diagonal vectors as each possible offset
		AB[i] = crypto.Zero(cps)
		for k := 0; k < n; k++ {
			// for row i, the kth diagonal contributes its value at column (i + k) % n
			offset := (i + k) % n
			maskedDiagonal := crypto.Mask(cps, diagonals[k][i], 0, false)
			// the maskedDiagonal has the diagonal value in its 0th index, so we need to rotate it over to the correct position before continuing
			// TODO maybe semantically we should do it in the other order? we can rotate over _then_ mask?
			maskedDiagonal = crypto.RotateRight(cps, maskedDiagonal, offset)

			// fmt.Printf("masked diagonal for i=%d, k=%d: ", i, k)
			// fmt.Print(crypto.DecryptFloatVector(cps, crypto.CipherVector{maskedDiagonal}, 2))
			// fmt.Print("\n")

			AB[i] = crypto.Add(cps, AB[i], maskedDiagonal)
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
