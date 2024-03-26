package main

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strconv"

	"github.com/BurntSushi/toml"
	"github.com/hhcho/sfgwas-private/gwas"
	"github.com/raulk/go-watchdog"

	mpc_core "github.com/hhcho/mpc-core"
	"math"
	"math/rand"
)

// Expects a party ID provided as an environment variable;
// e.g., run "PID=1 go run sfgwas.go"
var PID, PID_ERR = strconv.Atoi(os.Getenv("PID"))

// Default config path
var CONFIG_PATH = "config/"

func main() {
	// RunGWAS()
	RunMPCTest()
}

func InitProtocol(configPath string, mpcOnly bool) *gwas.ProtocolInfo {
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

	return gwas.InitializeGWASProtocol(config, PID, mpcOnly)
}

func RunGWAS() {
	if PID_ERR != nil {
		panic(PID_ERR)
	}

	// Initialize protocol
	prot := InitProtocol(CONFIG_PATH, false)

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

func RunMPCTest() {
	if PID_ERR != nil {
		panic(PID_ERR)
	}

	// Initialize protocol
	prot := InitProtocol(CONFIG_PATH, true)

	mpc := prot.GetMpc()[0]
	pid := mpc.GetPid()

	rtype := mpc_core.LElem256Zero
	fracBits := mpc.GetFracBits()

	scales := []float64{1e8, 1e4, 1e2, 1e1}
	N := uint(1e4)
	errorTol := 1e-4
  	
	for _, sc := range scales {
		x := make([]float64, N)
		answer := make([]float64, N)
		for i := range answer {
			x[i] = sc * rand.Float64()
			answer[i] = math.Sqrt(x[i])
		}

		var xRV mpc_core.RVec
		
		if pid == 1 { // Set shares to: x - mask

			xRV = mpc_core.FloatToRVec(rtype, x, fracBits)
		
			mpc.Network.Rand.SwitchPRG(2) // Shared PRG between 1 and 2
			mask := mpc.Network.Rand.RandVec(rtype, int(N))
			mpc.Network.Rand.RestorePRG()

			xRV.Sub(mask)

		} else if pid == 2 { // Set shares to: mask

			mpc.Network.Rand.SwitchPRG(1) // Shared PRG between 1 and 2
			mask := mpc.Network.Rand.RandVec(rtype, int(N))
			mpc.Network.Rand.RestorePRG()

			xRV = mask

		} else {
			xRV = mpc_core.InitRVec(rtype.Zero(), int(N))
		}

		s, _ := mpc.SqrtAndSqrtInverse(xRV, true)

		out := mpc.RevealSymVec(s).ToFloat(fracBits)

		if pid == 1 {
			any := false
			for j := range out {
				if (math.Abs(out[j] - answer[j]) / answer[j]) > errorTol {
					any = true
				}
			}

			if any {
				fmt.Printf("Sqrt test, scale %f: failed\n", sc)
			} else {
				fmt.Printf("Sqrt test, scale %f: success\n", sc)
			}
		}
	}

	prot.SyncAndTerminate(true)
}
