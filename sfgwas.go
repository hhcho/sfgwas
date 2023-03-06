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

	size := 4
	plaintextVector1 := make([]float64, size)
	plaintextVector2 := make([]float64, size)
	for i := 0; i < size; i++ {
		plaintextVector1[i] = float64(i*i)
		plaintextVector2[i] = float64(-i)
	}
	fmt.Println(plaintextVector1)
	fmt.Println(plaintextVector2)
	encryptedVector1, _ := crypto.EncryptFloatVector(selfParams, plaintextVector1)
	encryptedVector2, _ := crypto.EncryptFloatVector(selfParams, plaintextVector2)
	encSum := crypto.InnerSumAll(selfParams, encryptedVector1)
	decSum := crypto.DecryptFloat(selfParams, encSum)
	fmt.Println(decSum)

	encryptedResult := crypto.CAdd(selfParams, encryptedVector1, encryptedVector2) // CAdd standing for "Ciphertext Addition"
	// DecryptFloatVector both "Decrypts" (puts back into polynomial encoding) and "Decodes" (returns to message space)
	// note that, perhaps slightly confusingly, we call the polynomial encoding space the "plaintexts"
	plaintextResult := crypto.DecryptFloatVector(selfParams, encryptedResult, size)
	fmt.Println(plaintextResult)

	// multiply the matrices
	// 1 2 5      -2  0
	// 3 4 6       4 -1
	//             6  5
	// expected result:
	// 36 23
	// 46 26
	A := [][]float64{
		{1, 2, 5},
		{3, 4, 6},
	}
	// B := [][]float64{
	// 	{-2, 0},
	// 	{4, -1},
	// 	{6, 5},
	// }
	B := [][]float64{
		{1, -1},
		{1, 1},
		{1, 1},
	}
	// B := [][]float64 {
	// 	{0, -1, 5},
	// 	{-2, 4, 6},
	// }

	cipherA, _, _, _ := crypto.EncryptFloatMatrixRow(selfParams, A)
	m := len(B[0])

	cipherAB := cpMatrixMultiply(selfParams, cipherA, B)
	decAB := crypto.DecryptFloatMatrix(selfParams, cipherAB, m)
	fmt.Print(decAB)

	fmt.Print("sanity check")
	plainB, _, _, _ := crypto.EncodeFloatMatrixRow(selfParams, B)
	denseB := crypto.PlaintextToDense(selfParams, plainB, m)
	otherDecAB := gwas.CPMatMult1(selfParams, cipherA, denseB)
	fmt.Print(crypto.DecryptFloatMatrix(selfParams, otherDecAB, m))

	// TODO ask about bugs with InnerSum, InnerSumAll, etc. — are these params issues? does my current params setting lack rotation keys or something else needed? are those functions deprecated?
	// and maybe more broadly it's actually not immediately obvious how to perform a matrix multiplication precisely because of the difficulty of extracting the column vecgtor from a ciphermatrix (since each ciphervector has no internal structure that lets you pull arbitrary elements)

	// fmt.Println(len(A[0]))
	// fmt.Println(B[0])

	// AB := matrixMultiply(selfParams, A, B)
	// fmt.Println(AB)
	// cipherA, _, k1, _ := crypto.EncryptFloatMatrixRow(selfParams, A)
	// cipherB, _, m, _ := crypto.EncryptFloatMatrixRow(selfParams, B)
	// encryptedProd := cRowMatrixMultiply(selfParams, cipherA, cipherB, k1)
	// plaintextProd := crypto.DecryptFloatMatrix(selfParams, encryptedProd, m)
	// fmt.Println(plaintextProd)
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
