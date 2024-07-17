package gwas

import (
	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"github.com/ldsec/lattigo/v2/ckks"
	"go.dedis.ch/onet/v3/log"
	"math"
	"os"
	"strconv"
	"sync"
	"testing"
	"time"
)

// id of the node, by default we have two computing nodes and a coordinating party
var pid, _ = strconv.Atoi(os.Getenv("PID"))

func TestMicroBenchmark(t *testing.T) {
	// set MPC parameters and 3 nodes, default parameters
	servers := map[string]mpc.Server{}
	servers["party0"] = mpc.Server{
		IpAddr: "10.128.0.26", //"10.128.0.2", //"127.0.0.1",
		Ports:  map[string]string{"party1": "8010", "party2": "8020", "party3": "8030", "party4": "8040", "party5": "8050"},
	}
	servers["party1"] = mpc.Server{
		IpAddr: "10.128.0.5", //"10.128.0.5", //"127.0.0.1",
		Ports:  map[string]string{"party2": "8060", "party3": "8070", "party4": "8080", "party5": "8090"},
	}
	servers["party2"] = mpc.Server{
		IpAddr: "10.154.0.2", //"10.128.0.10", //"127.0.0.1",
		Ports:  map[string]string{"party3": "8100", "party4": "8110", "party5": "8120"},
	}
	/*servers["party3"] = mpc.Server{
		IpAddr: "10.128.0.15",
		Ports:  map[string]string{"party4": "8130", "party5": "8140"},
	}
	servers["party4"] = mpc.Server{
		IpAddr: "10.128.0.20",
		Ports:  map[string]string{"party5": "8150"},
	}
	servers["party5"] = mpc.Server{
		IpAddr: "10.128.0.21",
		Ports:  map[string]string{},
	}*/

	networks := mpc.ParallelNetworks(mpc.InitCommunication("0.0.0.0", servers, pid, 3, 1, ""))
	mpcEnv := mpc.InitParallelMPCEnv(networks, mpc_core.LElem256Zero, 60, 30)
	for thread := range mpcEnv {
		mpcEnv[thread].SetHubPid(1)
		mpcEnv[thread].SetBooleanShareFlag(true)
		mpcEnv[thread].SetDivSqrtMaxLen(10000)
	}

	// set MHE params, default parameters
	prec := uint(256)
	paramsMHE := ckks.DefaultParams[ckks.PN14QP438]
	for thread := range networks {
		networks[thread].SetMHEParams(paramsMHE)
	}

	cps := networks.CollectiveInit(paramsMHE, prec)

	// initiate variables
	encryptedVector := make(crypto.CipherVector, 1)
	encryptedVectorTwo := make(crypto.CipherVector, 1)
	plainVector := make(crypto.PlainVector, 1)

	vector := make([]float64, 8192)

	// Function to benchmark an operation 10 times
	benchmark := func(operation func() time.Duration, label string) {
		var totalDuration time.Duration
		networks.ResetNetworkLog()
		for i := 0; i < 10; i++ {
			totalDuration += operation()
			if i == 0 {
				networks.PrintNetworkLog()
				networks.ResetNetworkLog()
			}
		}
		avgDuration := totalDuration / 10
		log.LLvl1(pid, label, avgDuration)
	}

	if pid != 0 {
		//local operations: encryption and encoding (encoding only creates a plaintext)
		encryptedVectorTwo, _ = crypto.EncryptFloatVector(cps, vector)
		benchmark(func() time.Duration {
			timeBef := time.Now()
			encryptedVector, _ = crypto.EncryptFloatVector(cps, vector)
			return time.Since(timeBef)
		}, " Encryption time: ")
		benchmark(func() time.Duration {
			timeBef := time.Now()
			plainVector, _ = crypto.EncodeFloatVector(cps, vector)
			return time.Since(timeBef)
		}, " Encoding time: ")

		//local operations: addition, addition with plaintext, rotation
		benchmark(func() time.Duration {
			timeBef := time.Now()
			_ = crypto.CAdd(cps, encryptedVector, encryptedVector)
			return time.Since(timeBef)
		}, " Addition time: ")

		benchmark(func() time.Duration {
			timeBef := time.Now()
			_ = crypto.CPAdd(cps, encryptedVector, plainVector)
			return time.Since(timeBef)
		}, " Addition with plain time: ")

		benchmark(func() time.Duration {
			timeBef := time.Now()
			_ = crypto.RotateRight(cps, encryptedVector[0], 2)
			return time.Since(timeBef)
		}, " Rotation time: ")

		benchmark(func() time.Duration {
			timeBef := time.Now()
			_ = crypto.CPMult(cps, encryptedVector, plainVector)
			return time.Since(timeBef)
		}, " Multiplication with plain time: ")

		benchmark(func() time.Duration {
			timeBef := time.Now()
			_ = crypto.CMult(cps, encryptedVector, encryptedVector)
			return time.Since(timeBef)
		}, " Multiplication time: ")

		benchmark(func() time.Duration {
			timeBef := time.Now()
			_ = crypto.CMultScalar(cps, encryptedVector, encryptedVector[0])
			return time.Since(timeBef)
		}, " Multiplication with scalar time: ")

		benchmark(func() time.Duration {
			timeBef := time.Now()
			_ = crypto.InnerSumAll(cps, encryptedVector)
			return time.Since(timeBef)
		}, " Inner sum & Duplicate time: ")
	}

	mpcEnv[0].AssertSync()

	//distributed MHE operations

	// decryption
	benchmark(func() time.Duration {
		timeBef := time.Now()
		if pid != 0 {
			_ = networks[0].CollectiveDecrypt(cps, encryptedVector[0], 1)
		}
		timeFunc := time.Since(timeBef)
		mpcEnv[0].AssertSync()
		return timeFunc
	}, " Decryption time: ")

	benchmark(func() time.Duration {
		timeBef := time.Now()
		if pid != 0 {
			_ = networks[0].CollectiveBootstrapVec(cps, encryptedVector, 1)
		}
		timeFunc := time.Since(timeBef)
		mpcEnv[0].AssertSync()
		return timeFunc
	}, " Bootstrapping time: ")

	benchmark(func() time.Duration {
		timeBef := time.Now()
		if pid != 0 {
			_ = networks[0].AggregateCVec(cps, encryptedVectorTwo)
		}
		timeFunc := time.Since(timeBef)
		mpcEnv[0].AssertSync()
		return timeFunc
	}, " Aggregation time: ")

	benchmark(func() time.Duration {
		timeBef := time.Now()
		var mutex sync.Mutex
		if pid != 0 {
			_ = mpc.CSigmoidApprox(networks[0].GetPid(), networks[0], cps, encryptedVector,
				crypto.IntervalApprox{A: 0.1, B: 11.0, Degree: 15, Iter: 0, InverseNew: false}, &mutex)
		}
		timeFunc := time.Since(timeBef)
		mpcEnv[0].AssertSync()
		return timeFunc
	}, " Sigmoid time: ")

	var encryptedVectorSS, encryptedVectorSSInvSqrt mpc_core.RVec
	// switch to MPC (secret shares)
	benchmark(func() time.Duration {
		timeBef := time.Now()
		encryptedVectorSS = mpcEnv[0].CiphertextToSS(cps, mpcEnv[0].GetRType(),
			encryptedVector[0], 1, 8192)
		timeFunc := time.Since(timeBef)
		mpcEnv[0].AssertSync()
		return timeFunc
	}, " Switch to MPC time: ")
	mpcEnv[0].AssertSync()

	// inverse sqrt
	benchmark(func() time.Duration {
		timeBef := time.Now()
		_, encryptedVectorSSInvSqrt = mpcEnv[0].SqrtAndSqrtInverse(encryptedVectorSS, true)
		timeFunc := time.Since(timeBef)
		mpcEnv[0].AssertSync()
		return timeFunc
	}, " Inverse sqrt time: ")
	mpcEnv[0].AssertSync()

	// divide
	benchmark(func() time.Duration {
		timeBef := time.Now()
		_ = mpcEnv[0].Divide(encryptedVectorSS, encryptedVectorSS, true)
		timeFunc := time.Since(timeBef)
		mpcEnv[0].AssertSync()
		return timeFunc
	}, " Divide time: ")
	mpcEnv[0].AssertSync()

	// compare
	benchmark(func() time.Duration {
		timeBef := time.Now()
		_ = mpcEnv[0].IsPositive(encryptedVectorSS, true)
		timeFunc := time.Since(timeBef)
		mpcEnv[0].AssertSync()
		return timeFunc
	}, " Compare time: ")
	mpcEnv[0].AssertSync()

	rMat := mpc_core.InitRMat(encryptedVectorSS.Type().Zero(), 16, 16)

	// eigen decomposition
	//benchmark(func() time.Duration {
	networks.ResetNetworkLog()
	timeBef := time.Now()
	_, _ = mpcEnv[0].EigenDecomp(rMat)
	timeFunc := time.Since(timeBef)
	log.LLvl1("Eigendecomp: ", timeFunc)
	networks.PrintNetworkLog()
	mpcEnv[0].AssertSync()
	//	return timeFunc
	//}, " EigenDecomp time: ")
	//mpcEnv[0].AssertSync()

	// matrix multiplication
	networks.ResetNetworkLog()
	//benchmark(func() time.Duration {
	timeBef = time.Now()
	Vtr, Vtm := mpcEnv[0].BeaverPartitionMat(rMat)
	Lr, Lm := mpcEnv[0].BeaverPartitionMat(rMat)
	tempRes := mpcEnv[0].BeaverMultMat(Vtr, Vtm, Lr, Lm)
	_ = mpcEnv[0].BeaverReconstructMat(tempRes)
	timeFunc = time.Since(timeBef)
	log.LLvl1("Matrix mult time: ", timeFunc)
	networks.PrintNetworkLog()
	mpcEnv[0].AssertSync()
	//	return timeFunc
	//}, " Matrix mult time: ")
	mpcEnv[0].AssertSync()
	networks.ResetNetworkLog()

	//rMat = mpc_core.InitRMat(encryptedVectorSS.Type().Zero(), 16, 16)
	// matrix inverse sqrt SVD
	networks.ResetNetworkLog()
	//benchmark(func() time.Duration {
	timeBef = time.Now()
	invRMat := mpcEnv[0].MatrixInverseSqrtSVD(rMat)
	BT := mpcEnv[0].SSToCMat(cps, invRMat.Transpose())
	if pid != 0 {
		// scale back such that BTB ends up unscaled
		BT = crypto.CMultConstMat(cps, BT, math.Sqrt((1*2)/math.Sqrt(float64(16))), false)
		BT = crypto.CMatRescale(cps, BT)

		// Mask remaining values
		maskClear := make([]float64, 8192)
		for i := 0; i < 10; i++ {
			maskClear[i] = 1
		}
		maskEnc, _ := crypto.EncodeFloatVector(cps, maskClear)
		for i := range BT {
			BT[i] = crypto.CPMult(cps, BT[i], maskEnc)
		}

		_ = CMultMatInnerProd(cps, BT, crypto.CopyEncryptedMatrix(BT), 1)
	}
	timeFunc = time.Since(timeBef)
	log.LLvl1("MatrixInverseSqrtSVD time: ", timeFunc)
	networks.PrintNetworkLog()
	networks.ResetNetworkLog()
	mpcEnv[0].AssertSync()
	//return timeFunc
	//}, " MatrixInverseSqrtSVD time: ")
	mpcEnv[0].AssertSync()

	// switch back to MHE
	benchmark(func() time.Duration {
		timeBef := time.Now()
		_ = mpcEnv[0].SStoCiphertext(cps, encryptedVectorSSInvSqrt)
		timeFunc := time.Since(timeBef)
		mpcEnv[0].AssertSync()
		return timeFunc
	}, " Switch back to MHE time: ")
	mpcEnv[0].AssertSync()

}
