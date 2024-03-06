package mpc

import (
	"fmt"
	"gonum.org/v1/gonum/mat"
	"math"
	"os"
	"runtime"
	"sort"
	"strings"
	"sync"

	"github.com/aead/chacha20/chacha"
	"github.com/hhcho/frand"
	"github.com/hhcho/sfgwas/crypto"
	"github.com/ldsec/lattigo/v2/dckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"go.dedis.ch/onet/v3/log"

	"github.com/ldsec/lattigo/v2/ckks"
)

func (netObj ParallelNetworks) CollectiveInit(params *ckks.Parameters, prec uint) (cps *crypto.CryptoParams) {
	log.LLvl1("CollectiveInit started")

	dckksContext := dckks.NewContext(params)

	var kgen = ckks.NewKeyGenerator(params)

	var skShard *ckks.SecretKey
	if netObj[0].GetPid() == 0 {
		skShard = new(ckks.SecretKey)
		skShard.Value = dckksContext.RingQP.NewPoly()
	} else {
		skShard = kgen.GenSecretKey()

		skShard.Value.Zero()
		prng, err := utils.NewPRNG() // Use NewKeyedPRNG for debugging if deterministic behavior is desired
		if err != nil {
			panic(err)
		}
		ternarySamplerMontgomery := ring.NewTernarySampler(prng, dckksContext.RingQP, 1.0/3.0, true)
		skShard.Value = ternarySamplerMontgomery.ReadNew()
		dckksContext.RingQP.NTT(skShard.Value, skShard.Value)
	}

	// Globally shared random generator
	crpGen := make([]*ring.UniformSampler, len(netObj))
	for i := range crpGen {
		seed := make([]byte, chacha.KeySize)
		netObj[0].Rand.SwitchPRG(-1)
		netObj[0].Rand.RandRead(seed)
		netObj[0].Rand.RestorePRG()

		seedPrng := frand.NewCustom(seed, bufferSize, 20)

		crpGen[i] = ring.NewUniformSamplerWithBasePrng(seedPrng, dckksContext.RingQP)
	}

	log.LLvl1("PubKeyGen")
	var pk = netObj[0].CollectivePubKeyGen(params, skShard, crpGen[0])

	log.LLvl1("RelinKeyGen")
	var rlk = netObj[0].CollectiveRelinKeyGen(params, skShard, crpGen[0])

	nprocs := runtime.GOMAXPROCS(0)
	cps = crypto.NewCryptoParams(params, skShard, skShard, pk, rlk, prec, nprocs)

	smallDim := 20
	log.LLvl1("RotKeyGen: shifts <=", smallDim, "and powers of two up to", cps.GetSlots())
	if netObj[0].GetPid() > 0 {
		rotKs := netObj.CollectiveRotKeyGen(params, skShard, crpGen, crypto.GenerateRotKeys(cps.GetSlots(), smallDim, true))
		cps.RotKs = rotKs
		cps.SetEvaluators(cps.Params, rlk, cps.RotKs)
	}

	log.LLvl1("CollectiveInit finished")

	return
}

func (netObj *Network) CollectivePubKeyGen(parameters *ckks.Parameters, skShard *ckks.SecretKey, crpGen *ring.UniformSampler) (pk *ckks.PublicKey) {
	sk := &skShard.SecretKey

	ckgProtocol := dckks.NewCKGProtocol(parameters)

	pkShare := ckgProtocol.AllocateShares()

	crp := crpGen.ReadNew()
	ckgProtocol.GenShare(sk, crp, pkShare)

	pkAgg := netObj.AggregatePubKeyShares(pkShare)

	hubPid := netObj.GetHubPid()
	if netObj.GetPid() == 0 {
		pkAgg.Poly = netObj.ReceivePoly(hubPid)
	} else if netObj.GetPid() == hubPid {
		netObj.SendPoly(pkAgg.Poly, 0)
	}

	pk = ckks.NewPublicKey(parameters)
	ckgProtocol.GenPublicKey(pkAgg, crp, &pk.PublicKey)
	return
}

func (netObj *Network) CollectiveDecryptMat(cps *crypto.CryptoParams, cm crypto.CipherMatrix, sourcePid int) (pm crypto.PlainMatrix) {
	pid := netObj.GetPid()
	if pid == 0 {
		return
	}

	var tmp crypto.CipherMatrix
	var nr, nc int

	if sourcePid > 0 {
		if pid == sourcePid {
			nr, nc = len(cm), len(cm[0])
			for p := 1; p < netObj.GetNParty(); p++ {
				if p != sourcePid {
					netObj.SendInt(nr, p)
					netObj.SendInt(nc, p)
				}
			}
		} else {
			nr = netObj.ReceiveInt(sourcePid)
			nc = netObj.ReceiveInt(sourcePid)
			cm = make(crypto.CipherMatrix, nr)
			cm[0] = make(crypto.CipherVector, nc)
		}

		tmp = netObj.BroadcastCMat(cps, cm, sourcePid, nr, nc)
	} else {
		nr = len(cm)
		nc = len(cm[0])
		tmp = crypto.CopyEncryptedMatrix(cm)
	}

	tmp, level := crypto.FlattenLevels(cps, tmp)
	scale := tmp[0][0].Scale()

	parameters := cps.Params
	skShard := cps.Sk.Value

	zeroPoly := parameters.NewPolyQP()

	zeroPk := new(ckks.PublicKey)
	zeroPk.Value = [2]*ring.Poly{zeroPoly, zeroPoly}

	pcksProtocol := dckks.NewPCKSProtocol(parameters, 6.36)

	decShare := make([][]dckks.PCKSShare, nr)
	for i := range decShare {
		decShare[i] = make([]dckks.PCKSShare, nc)
		for j := range decShare[i] {
			decShare[i][j] = pcksProtocol.AllocateShares(level)
			pcksProtocol.GenShare(skShard, zeroPk, tmp[i][j], decShare[i][j])
		}
	}

	decAgg := netObj.AggregateDecryptSharesMat(decShare, level)

	pm = make(crypto.PlainMatrix, nr)
	for i := range pm {
		pm[i] = make(crypto.PlainVector, nc)
		for j := range pm[i] {
			ciphertextSwitched := ckks.NewCiphertext(parameters, 1, level, scale)
			pcksProtocol.KeySwitch(decAgg[i][j], tmp[i][j], ciphertextSwitched)
			pm[i][j] = ciphertextSwitched.Plaintext()
		}
	}

	return
}

func (netObj *Network) CollectiveDecryptVec(cps *crypto.CryptoParams, cv crypto.CipherVector, sourcePid int) (pv crypto.PlainVector) {
	if netObj.GetPid() == 0 {
		return
	}
	return netObj.CollectiveDecryptMat(cps, crypto.CipherMatrix{cv}, sourcePid)[0]
}

func (netObj *Network) CollectiveDecrypt(cps *crypto.CryptoParams, ct *ckks.Ciphertext, sourcePid int) (pt *ckks.Plaintext) {
	var tmp *ckks.Ciphertext

	// sourcePid broadcasts ct to other parties for collective decryption
	if netObj.GetPid() == sourcePid {
		for p := 1; p < netObj.GetNParty(); p++ {
			if p != sourcePid {
				netObj.SendCiphertext(ct, p)
			}
		}
		tmp = ct
	} else if netObj.GetPid() > 0 {
		tmp = netObj.ReceiveCiphertext(cps, sourcePid)
	} else { // pid == 0
		return
	}

	parameters := cps.Params
	skShard := cps.Sk.Value

	zeroPoly := parameters.NewPolyQP()

	zeroPk := new(ckks.PublicKey)
	zeroPk.Value = [2]*ring.Poly{zeroPoly, zeroPoly}

	pcksProtocol := dckks.NewPCKSProtocol(parameters, 6.36)

	decShare := pcksProtocol.AllocateShares(tmp.Level())
	pcksProtocol.GenShare(skShard, zeroPk, tmp, decShare)
	decAgg := netObj.AggregateDecryptShares(&decShare, tmp.Level())

	ciphertextSwitched := ckks.NewCiphertext(parameters, 1, tmp.Level(), tmp.Scale())
	pcksProtocol.KeySwitch(*decAgg, tmp, ciphertextSwitched)

	pt = ciphertextSwitched.Plaintext()

	return
}

func (netObj *Network) CollectiveBootstrap(cps *crypto.CryptoParams, ct *ckks.Ciphertext, sourcePid int) {
	// sourcePid broadcasts ct to other parties for collective decryption
	if netObj.GetPid() == 0 {
		return
	}

	// if sourcePid <= 0, assume cm is already shared across parties
	if sourcePid > 0 {
		if netObj.GetPid() == sourcePid {
			for p := 1; p < netObj.GetNParty(); p++ {
				if p != sourcePid {
					netObj.SendCiphertext(ct, p)
				}
			}
		} else {
			ct = netObj.ReceiveCiphertext(cps, sourcePid)
		}
	}

	parameters := cps.Params
	skShard := cps.Sk.Value
	crpGen := netObj.GetCRPGen()
	levelStart := ct.Level()

	refProtocol := dckks.NewRefreshProtocol(parameters)
	refShare1, refShare2 := refProtocol.AllocateShares(levelStart)

	crp := crpGen.ReadNew()

	refProtocol.GenShares(skShard, levelStart, netObj.GetNParty()-1, ct, parameters.Scale(), crp, refShare1, refShare2)

	refAgg1 := netObj.AggregateRefreshShare(refShare1, levelStart)
	refAgg2 := netObj.AggregateRefreshShare(refShare2, parameters.MaxLevel())

	refProtocol.Decrypt(ct, refAgg1)           // Masked decryption
	refProtocol.Recode(ct, parameters.Scale()) // Masked re-encoding
	refProtocol.Recrypt(ct, crp, refAgg2)      // Masked re-encryption

	return
}

func (netObj *Network) CollectiveBootstrapVec(cps *crypto.CryptoParams, cv crypto.CipherVector, sourcePid int) crypto.CipherVector {
	return netObj.CollectiveBootstrapMat(cps, crypto.CipherMatrix{cv}, sourcePid)[0]
}

func (netObj *Network) CollectiveBootstrapMat(cps *crypto.CryptoParams, cm crypto.CipherMatrix, sourcePid int) crypto.CipherMatrix {
	if netObj.GetPid() == 0 {
		return cm
	}

	// if sourcePid <= 0, assume cm is already shared across parties
	if sourcePid > 0 {
		// sourcePid broadcasts ct to other parties for collective decryption
		if netObj.GetPid() == sourcePid {
			for p := 1; p < netObj.GetNParty(); p++ {
				if p != sourcePid {
					netObj.SendInt(len(cm), p)
					netObj.SendInt(len(cm[0]), p)
					netObj.SendCipherMatrix(cm, p)
				}
			}
		} else {
			ncols := netObj.ReceiveInt(sourcePid)
			nrows := netObj.ReceiveInt(sourcePid)
			cm = netObj.ReceiveCipherMatrix(cps, ncols, nrows, sourcePid)
		}
	}

	parameters := cps.Params
	skShard := cps.Sk.Value
	crpGen := netObj.GetCRPGen()

	cm, levelStart := crypto.FlattenLevels(cps, cm)
	//log.LLvl1(time.Now().Format(time.RFC3339), "Bootstrap: dimensions", len(cm), "x", len(cm[0]), "input level", levelStart)

	refProtocol := dckks.NewRefreshProtocol(parameters)

	refSharesDecrypt := make([][]ring.Poly, len(cm))
	crps := make([][]ring.Poly, len(cm))
	refSharesRecrypt := make([][]ring.Poly, len(cm))

	for i := range cm {
		refSharesDecrypt[i] = make([]ring.Poly, len(cm[0]))
		crps[i] = make([]ring.Poly, len(cm[0]))
		refSharesRecrypt[i] = make([]ring.Poly, len(cm[0]))

		for j := range cm[i] {
			refShare1, refShare2 := refProtocol.AllocateShares(levelStart)

			refSharesDecrypt[i][j] = *refShare1
			refSharesRecrypt[i][j] = *refShare2
			crps[i][j] = *crpGen.ReadNew()

			refProtocol.GenShares(skShard, levelStart, netObj.GetNParty()-1, cm[i][j], parameters.Scale(),
				&(crps[i][j]), &(refSharesDecrypt[i][j]), &(refSharesRecrypt[i][j]))

		}
	}

	refAgg1 := netObj.AggregateRefreshShareMat(refSharesDecrypt, levelStart)
	refAgg2 := netObj.AggregateRefreshShareMat(refSharesRecrypt, parameters.MaxLevel())

	for i := range cm {
		for j := range cm[i] {
			//no communication
			refProtocol.Decrypt(cm[i][j], &refAgg1[i][j])              // Masked decryption
			refProtocol.Recode(cm[i][j], parameters.Scale())           // Masked re-encoding
			refProtocol.Recrypt(cm[i][j], &crps[i][j], &refAgg2[i][j]) // Masked re-encryption

			// Fix discrepancy in number of moduli
			if len(cm[i][j].Value()[0].Coeffs) < len(cm[i][j].Value()[1].Coeffs) {
				poly := ring.NewPoly(len(cm[i][j].Value()[0].Coeffs[0]), len(cm[i][j].Value()[0].Coeffs))
				for pi := range poly.Coeffs {
					for pj := range poly.Coeffs[0] {
						poly.Coeffs[pi][pj] = cm[i][j].Value()[1].Coeffs[pi][pj]
					}
				}
				cm[i][j].Value()[1] = poly
			}
		}
	}

	//log.LLvl1(time.Now().Format(time.RFC3339), "Bootstrap: output level", cm[0][0].Level())

	return cm

}

// BootstrapMatAll: collective bootstrap for all parties (except 0)
func (netObj *Network) BootstrapMatAll(cps *crypto.CryptoParams, cm crypto.CipherMatrix) crypto.CipherMatrix {

	tmp := make(crypto.CipherMatrix, len(cm))

	//TODO: optimize to run simultaneously
	for sourcePid := 1; sourcePid < netObj.GetNParty(); sourcePid++ {
		if netObj.GetPid() == sourcePid {
			cm = netObj.CollectiveBootstrapMat(cps, cm, sourcePid)
		} else {
			netObj.CollectiveBootstrapMat(cps, tmp, sourcePid)
		}
	}

	return cm
}

func (netObj *Network) BootstrapVecAll(cps *crypto.CryptoParams, cv crypto.CipherVector) crypto.CipherVector {
	tmp := make(crypto.CipherVector, len(cv))

	for sourcePid := 1; sourcePid < netObj.GetNParty(); sourcePid++ {
		if netObj.GetPid() == sourcePid {
			cv = netObj.CollectiveBootstrapVec(cps, cv, sourcePid)
		} else {
			netObj.CollectiveBootstrapVec(cps, tmp, sourcePid)
		}
	}

	return cv
}

func (netObj ParallelNetworks) CollectiveRotKeyGen(parameters *ckks.Parameters, skShard *ckks.SecretKey,
	crpGen []*ring.UniformSampler, rotTypes []crypto.RotationType) (rotKeys *ckks.RotationKeySet) {

	slots := parameters.Slots()

	sk := &skShard.SecretKey

	shiftMap := make(map[int]bool)
	for _, rotType := range rotTypes {

		var shift int
		if rotType.Side == crypto.SideRight {
			shift = slots - rotType.Value
		} else {
			shift = rotType.Value
		}

		shiftMap[shift] = true
	}

	gElems := make([]uint64, len(shiftMap))
	i := 0
	for k := range shiftMap {
		gElems[i] = parameters.GaloisElementForColumnRotationBy(k)
		i++
	}
	// Add rotation for taking the complex conjugate
	gElems = append(gElems, parameters.GaloisElementForRowRotation())

	// Need to sortInt otherwise different parties might have different ordering
	sort.Slice(gElems, func(i, j int) bool { return gElems[i] < gElems[j] })

	rotKeys = ckks.NewRotationKeySet(parameters, gElems)
	//
	//for ind, galEl := range gElems {
	//	rtgProtocol := dckks.NewRotKGProtocol(parameters)
	//	rtgShare := rtgProtocol.AllocateShares()
	//
	//	rtgCrp := make([]*ring.Poly, parameters.Beta())
	//	for i := 0; i < parameters.Beta(); i++ {
	//		rtgCrp[i] = crpGen[0].ReadNew()
	//	}
	//
	//	rtgProtocol.GenShare(sk, galEl, rtgCrp, rtgShare)
	//
	//	rtgAgg := netObj[0].AggregateRotKeyShare(rtgShare)
	//
	//	rtgProtocol.GenRotationKey(rtgAgg, rtgCrp, rotKeys.Keys[galEl])
	//	fmt.Println("Generate RotKey ", ind+1, "/", len(gElems), ", Galois element", galEl)
	//}

	/* Parallel version */
	nproc := len(netObj)
	jobChannels := make([]chan uint64, nproc)
	for i := range jobChannels {
		jobChannels[i] = make(chan uint64, 32)
	}

	// Dispatcher
	go func() {
		for ind, galEl := range gElems {
			jobChannels[ind%nproc] <- galEl
			fmt.Println("Generate RotKey ", ind+1, "/", len(gElems), ", Galois element", galEl)
		}
		for _, c := range jobChannels {
			close(c)
		}
	}()

	// Workers
	var wg sync.WaitGroup
	for thread := 0; thread < nproc; thread++ {
		wg.Add(1)
		go func(thread int, net *Network, crpGen *ring.UniformSampler) {
			defer wg.Done()
			for galEl := range jobChannels[thread] {
				rtgProtocol := dckks.NewRotKGProtocol(parameters)
				rtgShare := rtgProtocol.AllocateShares()

				rtgCrp := make([]*ring.Poly, parameters.Beta())
				for i := 0; i < parameters.Beta(); i++ {
					rtgCrp[i] = crpGen.ReadNew()
				}

				rtgProtocol.GenShare(sk, galEl, rtgCrp, rtgShare)

				rtgAgg := net.AggregateRotKeyShare(rtgShare)

				rtgProtocol.GenRotationKey(rtgAgg, rtgCrp, rotKeys.Keys[galEl])
			}
		}(thread, netObj[thread], crpGen[thread])
	}
	wg.Wait()

	return
}

func (netObj *Network) CollectiveRelinKeyGen(params *ckks.Parameters, skShard *ckks.SecretKey, crpGen *ring.UniformSampler) (evk *ckks.RelinearizationKey) {
	sk := &skShard.SecretKey

	prot := dckks.NewRKGProtocol(params)
	ephSk, share1, share2 := prot.AllocateShares()

	crp := make([]*ring.Poly, params.Beta())
	for i := 0; i < params.Beta(); i++ {
		crp[i] = crpGen.ReadNew()
	}

	evk = ckks.NewRelinearizationKey(params)

	if netObj.GetPid() > 0 {
		prot.GenShareRoundOne(sk, crp, ephSk, share1)
		outRound1 := netObj.AggregateRelinKeyShare(share1, true)

		prot.GenShareRoundTwo(ephSk, sk, outRound1, crp, share2)
		outRound2 := netObj.AggregateRelinKeyShare(share2, true)

		prot.GenRelinearizationKey(outRound1, outRound2, &evk.RelinearizationKey)
	}

	return
}

func (netObj *Network) BroadcastCMat(cps *crypto.CryptoParams, cm crypto.CipherMatrix, sourcePid int, numCtxRow, numCtxCol int) crypto.CipherMatrix {
	if netObj.GetPid() == sourcePid {
		if len(cm) != numCtxRow || len(cm[0]) != numCtxCol {
			panic("BroadcastCVec: dimensions of cm do not match numCtxRow or numCtxCol")
		}

		for p := 1; p < netObj.GetNParty(); p++ {
			if p != sourcePid {
				netObj.SendCipherMatrix(cm, p)
			}
		}
		cm = crypto.CopyEncryptedMatrix(cm)
	} else if netObj.GetPid() > 0 {
		cm = netObj.ReceiveCipherMatrix(cps, numCtxRow, numCtxCol, sourcePid)
	}
	return cm

}
func (netObj *Network) BroadcastCVec(cps *crypto.CryptoParams, cv crypto.CipherVector, sourcePid int, numCtx int) crypto.CipherVector {
	if netObj.GetPid() == sourcePid {
		if len(cv) != numCtx {
			panic("BroadcastCVec: len(cv) does not match numCtx")
		}

		for p := 1; p < netObj.GetNParty(); p++ {
			if p != sourcePid {
				netObj.SendCipherVector(cv, p)
			}
		}
		cv = crypto.CopyEncryptedVector(cv)
	} else if netObj.GetPid() > 0 {
		cv = netObj.ReceiveCipherVector(cps, numCtx, sourcePid)
	}
	return cv
}

func (netObj *Network) BroadcastCiphertext(cps *crypto.CryptoParams, ct *ckks.Ciphertext, sourcePid int) *ckks.Ciphertext {
	if netObj.GetPid() == sourcePid {
		for p := 1; p < netObj.GetNParty(); p++ {
			if p != sourcePid {
				netObj.SendCiphertext(ct, p)
			}
		}
		ct = ct.CopyNew().Ciphertext()
	} else if netObj.GetPid() > 0 {
		ct = netObj.ReceiveCiphertext(cps, sourcePid)
	}
	return ct
}

func SaveMatrixToFileWithPrint(cps *crypto.CryptoParams, mpcObj *MPC, cm crypto.CipherMatrix, nElemCol int, sourcePid int, filename string, print bool) {
	SaveMatrixToFileWithPrintIndex(cps, mpcObj, cm, nElemCol, sourcePid, filename, print, 0, 10)
}

// SaveMatrixToFileWithPrint saves a matrix to a file and prints the first nElemCol elements of each row
func SaveMatrixToFileWithPrintIndex(cps *crypto.CryptoParams, mpcObj *MPC, cm crypto.CipherMatrix, nElemCol int, sourcePid int, filename string, print bool, firstIndex, lastIndex int) {
	log.LLvl1("Saving matrix to file", filename, "from pid", sourcePid, "with", len(cm), "rows and", nElemCol, "columns")
	pid := mpcObj.GetPid()
	if pid == 0 {
		return
	}

	pm := mpcObj.Network.CollectiveDecryptMat(cps, cm, sourcePid)

	if pid == sourcePid || sourcePid < 0 {

		M := mat.NewDense(len(cm), nElemCol, nil)
		for i := range pm {
			decodedRow := crypto.DecodeFloatVector(cps, pm[i])[:nElemCol]
			if print {
				log.LLvl1(filename, " : ", decodedRow[firstIndex:int(crypto.Min(nElemCol, lastIndex))])
			}

			M.SetRow(i, decodedRow)
		}
		log.LLvl1("Matrix decoded")

		f, err := os.Create(filename)

		if err != nil {
			panic(err)
		}

		defer f.Close()

		rows, cols := M.Dims()

		for row := 0; row < rows; row++ {
			line := make([]string, cols)
			for col := 0; col < cols; col++ {
				line[col] = fmt.Sprintf("%.6e", M.At(row, col))
			}

			f.WriteString(strings.Join(line, ",") + "\n")
		}

		f.Sync()

		fmt.Println("Saved data to", filename)

	}

}

func CSigmoidApprox(sourcePid int, net *Network, cryptoParams *crypto.CryptoParams, ctIn crypto.CipherVector,
	intv crypto.IntervalApprox, mutex *sync.Mutex) crypto.CipherVector {
	res := make(crypto.CipherVector, len(ctIn))
	for i := range ctIn {
		res[i] = SigmoidApprox(sourcePid, net, cryptoParams, ctIn[i], intv, mutex)
	}
	return res
}

func SigmoidApprox(sourcePid int, net *Network, cryptoParams *crypto.CryptoParams, ctIn *ckks.Ciphertext,
	intv crypto.IntervalApprox, mutex *sync.Mutex) *ckks.Ciphertext {
	var y *ckks.Ciphertext

	if intv.Degree == 0 {
		ctDecrypt := net.CollectiveDecrypt(cryptoParams, ctIn, sourcePid)
		cdDecode := crypto.DecodeFloatVector(cryptoParams, crypto.PlainVector{ctDecrypt})
		cdOut := make([]float64, len(cdDecode))
		log.LLvl1("before approx result ", cdDecode[:300])
		for i := 0; i < len(cdDecode); i++ {
			cdOut[i] = 1.0 / math.Sqrt(cdDecode[i])
		}
		log.LLvl1("approx result ", cdOut[:300])
		cdOutEncrypted, _ := crypto.EncryptFloatVector(cryptoParams, cdOut)
		return cdOutEncrypted[0]
	} else {

		cheby := ckks.Approximate(Sigmoid, complex(intv.A, 0), complex(intv.B, 0), intv.Degree)
		// We evaluate the interpolated Chebyshev interpolant on y
		// Change of variable
		a := cheby.A()
		b := cheby.B()

		if ctIn.Level() < int(math.Ceil(math.Log2(float64(intv.Degree+1)))+1)+2 {
			mutex.Lock()
			ctIn = net.BootstrapVecAll(cryptoParams, crypto.CipherVector{ctIn})[0]
			mutex.Unlock()
		}

		err := cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			y = eval.MultByConstNew(ctIn.CopyNew().Ciphertext(), 2/(b-a))
			if err := eval.Rescale(y, cryptoParams.Params.Scale(), y); err != nil {
				panic(err)
			}
			eval.AddConst(y, (-a-b)/(b-a), y)
			return nil
		})
		if err != nil {
			log.Fatal(err)
		}

		ctInCopy := y.CopyNew().Ciphertext()

		err = cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			var err error
			y, err = eval.EvaluateCheby(ctInCopy, cheby, ctInCopy.Scale())
			if err != nil {
				log.Fatal(err)
			}
			return err
		})
		if err != nil {
			log.Fatal(err)
		}
	}
	return y
}

func Sigmoid(x complex128) complex128 {
	return complex(1.0/(1+math.Exp(-real(x))), 0)
}
