package mpc

import (
	mpc_core "github.com/hhcho/mpc-core"
)

func (mpcObj *MPC) BeaverPartition(a mpc_core.RElem) (mpc_core.RElem, mpc_core.RElem) {
	ar, am := mpcObj.BeaverPartitionMat(mpc_core.RMat{mpc_core.RVec{a}})
	return ar[0][0], am[0][0]
}

func (mpcObj *MPC) BeaverPartitionVec(a mpc_core.RVec) (mpc_core.RVec, mpc_core.RVec) {
	ar, am := mpcObj.BeaverPartitionMat(mpc_core.RMat{a})
	return ar[0], am[0]
}

// Returns ar, am
func (mpcObj *MPC) BeaverPartitionMat(a mpc_core.RMat) (mpc_core.RMat, mpc_core.RMat) {
	pid := mpcObj.Network.pid
	nrows, ncols := a.Dims()
	rtype := a.Type()

	var mask mpc_core.RMat

	if pid == 0 {
		am := mpc_core.InitRMat(rtype.Zero(), nrows, ncols)

		for p := 1; p < mpcObj.Network.NumParties; p++ {
			mpcObj.Network.Rand.SwitchPRG(p)
			mask = mpcObj.Network.Rand.RandMat(rtype, nrows, ncols)
			am.Add(mask)
			mpcObj.Network.Rand.RestorePRG()
		}

		return mpc_core.InitRMat(rtype.Zero(), nrows, ncols), am
	}

	// Retrieve random shares from pid = 0
	mpcObj.Network.Rand.SwitchPRG(0)
	mask = mpcObj.Network.Rand.RandMat(rtype, nrows, ncols)
	mpcObj.Network.Rand.RestorePRG()

	ar := a.Copy()
	ar.Sub(mask)
	ar = mpcObj.RevealSymMat(ar)
	return ar, mask
}

func (mpcObj *MPC) BeaverReconstruct(a mpc_core.RElem) mpc_core.RElem {
	return mpcObj.BeaverReconstructMat(mpc_core.RMat{mpc_core.RVec{a}})[0][0]
}

func (mpcObj *MPC) BeaverReconstructVec(a mpc_core.RVec) mpc_core.RVec {
	return mpcObj.BeaverReconstructMat(mpc_core.RMat{a})[0]
}

func (mpcObj *MPC) BeaverReconstructMat(a mpc_core.RMat) mpc_core.RMat {
	pid := mpcObj.Network.pid

	rtype := a.Type()
	nr, nc := a.Dims()

	last := mpcObj.Network.NumParties - 1

	if pid == 0 {

		mask := a.Copy()
		for to := 1; to < mpcObj.Network.NumParties-1; to++ {
			mpcObj.Network.Rand.SwitchPRG(to)
			share := mpcObj.Network.Rand.RandMat(rtype, nr, nc)
			mpcObj.Network.Rand.RestorePRG()
			mask.Sub(share)
		}
		mpcObj.Network.SendRData(mask, last)
		return mask //return a

	}

	var mask mpc_core.RMat
	if pid == last {
		mask = mpcObj.Network.ReceiveRMat(rtype, nr, nc, 0)
	} else {
		mpcObj.Network.Rand.SwitchPRG(0)
		mask = mpcObj.Network.Rand.RandMat(rtype, nr, nc)
		mpcObj.Network.Rand.RestorePRG()
	}

	ar := a.Copy()
	ar.Add(mask)

	return ar
}

func (mpcObj *MPC) BeaverMult(ar, am, br, bm mpc_core.RElem) mpc_core.RElem {
	// this function is not used by sf-relate
	pid := mpcObj.Network.pid
	if pid == 0 {
		return am.Mul(bm)
	}

	out := ar.Mul(bm)
	out = out.Add(br.Mul(am))
	if pid == 1 {
		out = out.Add(ar.Mul(br))
	}
	return out
}

func (mpcObj *MPC) BeaverMultElemVec(ar, am, br, bm mpc_core.RVec) mpc_core.RVec {
	return mpcObj.BeaverMultElemMat(mpc_core.RMat{ar}, mpc_core.RMat{am}, mpc_core.RMat{br}, mpc_core.RMat{bm}, mpcObj.Network.pid == 1)[0]
}

func (mpcObj *MPC) BeaverMultElemMat(ar, am, br, bm mpc_core.RMat, sending bool) mpc_core.RMat {
	pid := mpcObj.Network.pid
	nr, nc := am.Dims()

	if pid == 0 {
		out := am.Copy()
		out.MulElem(bm)
		return out
	}

	out := mpc_core.InitRMat(am.Type().Zero(), nr, nc)
	for i := 0; i < nr; i++ {
		for j := 0; j < nc; j++ {
			out[i][j] = out[i][j].Add(ar[i][j].Mul(bm[i][j]))
			out[i][j] = out[i][j].Add(br[i][j].Mul(am[i][j]))
			// ensure only one party is doing this
			if sending {
				out[i][j] = out[i][j].Add(ar[i][j].Mul(br[i][j]))
			}
		}
	}
	return out
}

func (mpcObj *MPC) BeaverMultMat(ar, am, br, bm mpc_core.RMat) mpc_core.RMat {
	pid := mpcObj.Network.pid
	if pid == 0 {
		return mpc_core.RMultMat(am, bm)
	}

	out := mpc_core.RMultMat(ar, bm)
	out.Add(mpc_core.RMultMat(am, br))
	if pid == 1 {
		out.Add(mpc_core.RMultMat(ar, br))
	}
	return out
}
