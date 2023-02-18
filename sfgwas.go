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

func playground() {
	precision := uint(8) // 8 and 32 seem to work well, while 16 for some reason seems to accumulate fairly large error (10^-3 ish)
	params := crypto.NewCryptoParamsForNetwork(ckks.DefaultParams[ckks.PN13QP218], 1, precision)
	selfParams := params[0]

	size := 4
	plaintextVector1 := make([]float64, size)
	plaintextVector2 := make([]float64, size)
	for i := 0; i < size; i++ {
		plaintextVector1[i] = float64(i*i)
		plaintextVector2[i] = float64(-i)
	}
	encryptedVector1, _ := crypto.EncryptFloatVector(selfParams, plaintextVector1)
	encryptedVector2, _ := crypto.EncryptFloatVector(selfParams, plaintextVector2)

	encryptedResult := crypto.CAdd(selfParams, encryptedVector1, encryptedVector2) // CAdd standing for "Ciphertext Addition"
	// DecryptFloatVector both "Decrypts" (puts back into polynomial encoding) and "Decodes" (returns to message space)
	// note that, perhaps slightly confusingly, we call the polynomial encoding space the "plaintexts"
	plaintextResult := crypto.DecryptFloatVector(selfParams, encryptedResult, size)
	fmt.Print(plaintextResult)
}
