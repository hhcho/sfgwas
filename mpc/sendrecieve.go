package mpc

import (
	"encoding/binary"
	"fmt"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas-private/crypto"
	"github.com/ldsec/lattigo/v2/ring"

	"github.com/ldsec/lattigo/v2/ckks"
)

// TODO signedness lost
func (netObj *Network) SendInt(val, to int) {
	conn := netObj.conns[to]

	buf := make([]byte, 8)
	binary.LittleEndian.PutUint64(buf, uint64(val))
	WriteFull(&conn, buf)

	netObj.UpdateSenderLog(to, 8)

}

func (netObj *Network) SendIntVector(v []uint64, to int) {
	conn := netObj.conns[to]

	n := uint64(len(v))
	bytes := make([]byte, 8*n)
	for i := range v {
		binary.LittleEndian.PutUint64(bytes[(i*8):((i*8)+8)], v[i])
	}

	WriteFull(&conn, bytes)

	netObj.UpdateSenderLog(to, len(bytes))
}

//SendCiphertext sends ciphertext over a connection
func (netObj *Network) SendCiphertext(ct *ckks.Ciphertext, to int) {
	conn := netObj.conns[to]

	bytes, e := ct.MarshalBinary() //convert to bytes
	checkError(e)

	sbuf := make([]byte, 8)
	binary.LittleEndian.PutUint64(sbuf, uint64(len(bytes)))

	WriteFull(&conn, sbuf)
	WriteFull(&conn, bytes)

	netObj.UpdateSenderLog(to, 8+len(bytes))

	return
}

//SendCipherVector sends ciphervector over a connection
func (netObj *Network) SendCipherVector(v crypto.CipherVector, to int) {
	conn := netObj.conns[to]

	sbytes, cvbytes := MarshalCV(v)

	sbuf := make([]byte, 8) //TODO: see if 4 bytes is enough
	cvbuf := make([]byte, 8)
	binary.LittleEndian.PutUint64(sbuf, uint64(len(sbytes)))
	binary.LittleEndian.PutUint64(cvbuf, uint64(len(cvbytes)))

	WriteFull(&conn, sbuf)
	WriteFull(&conn, sbytes)
	WriteFull(&conn, cvbuf)
	WriteFull(&conn, cvbytes)

	netObj.UpdateSenderLog(to, len(sbytes)+len(cvbytes))

}

//SendCipherMatrix sends ciphermatrix over a connection
func (netObj *Network) SendCipherMatrix(m crypto.CipherMatrix, to int) {
	conn := netObj.conns[to]

	sbytes, cmbytes := crypto.MarshalCM(m)

	sbuf := make([]byte, 8) //TODO: see if 4 bytes is enough
	cmbuf := make([]byte, 8)
	binary.LittleEndian.PutUint64(sbuf, uint64(len(sbytes)))
	binary.LittleEndian.PutUint64(cmbuf, uint64(len(cmbytes)))

	WriteFull(&conn, sbuf)
	WriteFull(&conn, sbytes)
	WriteFull(&conn, cmbuf)
	WriteFull(&conn, cmbytes)

	netObj.UpdateSenderLog(to, len(sbytes)+len(cmbytes))
}

//SendAllCryptoParams sends cryptoParams to all parties (receives nothing)
func (netObj *Network) SendAllCryptoParams(cps *crypto.CryptoParams) {

	buf := make([]byte, 8)
	bytes, _ := cps.MarshalBinary()
	binary.LittleEndian.PutUint64(buf, uint64(len(bytes)))

	for i := 0; i < netObj.NumParties; i++ {
		if i != netObj.pid {

			conn := netObj.conns[i]
			WriteFull(&conn, buf)
			WriteFull(&conn, bytes)
		}
	}
}

//SendAllCiphertext sends ciphertext over to all parties
func (netObj *Network) SendAllCiphertext(ct *ckks.Ciphertext, includeZero bool) {
	for i := 0; i < netObj.NumParties; i++ {
		if i == 0 && !includeZero {
			continue
		}
		if i != netObj.pid {
			netObj.SendCiphertext(ct, i)
		}
	}
}

//SendAllCipherVector sends ciphevector over to all parties
func (netObj *Network) SendAllCipherVector(cv crypto.CipherVector, includeZero bool) {
	for i := 0; i < netObj.NumParties; i++ {
		if i == 0 && !includeZero {
			continue
		}
		if i != netObj.pid {
			netObj.SendCipherVector(cv, i)
		}
	}
}

//SendAllCipherMatrix sends ciphermatrix over to all parties
func (netObj *Network) SendAllCipherMatrix(cm crypto.CipherMatrix, includeZero bool) {
	for i := 0; i < netObj.NumParties; i++ {
		if i == 0 && !includeZero {
			continue
		}
		if i != netObj.pid {
			netObj.SendCipherMatrix(cm, i)
		}
	}
}

func (netObj *Network) SendPoly(poly *ring.Poly, pid int) {
	conn := netObj.conns[pid]
	bytes, _ := poly.MarshalBinary() //convert to bytes

	buf := make([]byte, 4)
	binary.LittleEndian.PutUint32(buf, uint32(len(bytes)))

	WriteFull(&conn, buf)
	WriteFull(&conn, bytes)

	netObj.UpdateSenderLog(pid, len(buf)+len(bytes))
}

func (netObj *Network) SendPolyMat(poly [][]ring.Poly, pid int) {
	conn := netObj.conns[pid]
	//marshal poly
	sizes, bytes := MarshalPolyMat(poly)

	//send over sizes
	sbuf := make([]byte, 8)

	binary.LittleEndian.PutUint64(sbuf, uint64(len(sizes)))
	WriteFull(&conn, sbuf)
	WriteFull(&conn, sizes)

	binary.LittleEndian.PutUint64(sbuf, uint64(len(bytes)))
	WriteFull(&conn, sbuf)
	WriteFull(&conn, bytes)

	netObj.UpdateSenderLog(pid, 2*len(sbuf)+len(sizes)+len(bytes))

}

//SendRData sends ring/field elements to party p
func (netObj *Network) SendRData(data interface{}, p int) {
	conn := netObj.conns[p]
	bytes := MarshalRData(data) //convert to bytes

	var buf []byte

	switch data.(type) {
	case mpc_core.RElem:
	default:
		buf = make([]byte, 4)
		binary.LittleEndian.PutUint32(buf, uint32(len(bytes)))
		WriteFull(&conn, buf)
	}

	WriteFull(&conn, bytes)

	netObj.UpdateSenderLog(p, len(buf)+len(bytes))
}

/*ALL parties need to share information*/

func (netObj *Network) ReceiveInt(from int) int {
	conn := netObj.conns[from]

	buf := make([]byte, 8)
	ReadFull(&conn, buf)
	val := int(binary.LittleEndian.Uint64(buf))

	netObj.UpdateReceiverLog(from, len(buf))

	return val
}

//ReceiveCipherVector reads bytes sent over conn, and unmarshals it as a ciphervector
func (netObj *Network) ReceiveIntVector(nElem, from int) []uint64 {
	conn := netObj.conns[from]

	data := make([]byte, nElem*8)
	ReadFull(&conn, data)

	out := make([]uint64, nElem)
	for i := range out {
		out[i] = binary.LittleEndian.Uint64(data[(i * 8):((i * 8) + 8)])
	}

	netObj.UpdateReceiverLog(from, len(data))

	return out
}

// // ReceiveCiphertext reads bytes sent over conn, and unmarshals it as a ciphertext
func (netObj *Network) ReceiveCiphertext(cryptoParams *crypto.CryptoParams, from int) *ckks.Ciphertext {
	conn := netObj.conns[from]

	buf := make([]byte, 8)
	ReadFull(&conn, buf)
	sbyteSize := binary.LittleEndian.Uint64(buf)
	data := make([]byte, sbyteSize)
	ReadFull(&conn, data)

	ct := ckks.NewCiphertext(cryptoParams.Params, 1, cryptoParams.Params.MaxLevel(), cryptoParams.Params.Scale())
	e := ct.UnmarshalBinary(data)
	checkError(e)

	netObj.UpdateReceiverLog(from, len(buf)+len(data))
	return ct
}

//ReceiveCipherVector reads bytes sent over conn, and unmarshals it as a ciphervector
func (netObj *Network) ReceiveCipherVector(cryptoParams *crypto.CryptoParams, ncts, from int) crypto.CipherVector {
	conn := netObj.conns[from]

	sbuf := make([]byte, 8)
	ReadFull(&conn, sbuf)

	sbyteSize := binary.LittleEndian.Uint64(sbuf)
	sdata := make([]byte, sbyteSize)
	ReadFull(&conn, sdata)

	cmbuf := make([]byte, 8)
	ReadFull(&conn, cmbuf)

	cbyteSize := binary.LittleEndian.Uint64(cmbuf)
	cdata := make([]byte, cbyteSize)
	ReadFull(&conn, cdata)

	netObj.UpdateReceiverLog(from, len(sbuf)+len(cmbuf)+len(cdata)+len(sdata))

	return UnmarshalCV(cryptoParams, ncts, sdata, cdata)
}

//ReceiveCipherMatrix reads bytes sent over conn, and unmarshals it as a ciphermatrix
func (netObj *Network) ReceiveCipherMatrix(cryptoParams *crypto.CryptoParams, nv, nct, from int) crypto.CipherMatrix {
	conn := netObj.conns[from]

	sbuf := make([]byte, 8)
	ReadFull(&conn, sbuf)
	sbyteSize := binary.LittleEndian.Uint64(sbuf)
	sdata := make([]byte, sbyteSize)
	ReadFull(&conn, sdata)

	cmbuf := make([]byte, 8)
	ReadFull(&conn, cmbuf)
	cbyteSize := binary.LittleEndian.Uint64(cmbuf)
	cdata := make([]byte, cbyteSize)
	ReadFull(&conn, cdata)

	netObj.UpdateReceiverLog(from, len(sbuf)+len(sdata)+len(cmbuf)+len(cdata))
	return crypto.UnmarshalCM(cryptoParams, nv, nct, sdata, cdata)
}

func (netObj *Network) ReceivePoly(pid int) *ring.Poly {
	conn := netObj.conns[pid]

	buf := make([]byte, 4)
	ReadFull(&conn, buf)

	byteSize := binary.LittleEndian.Uint32(buf)
	data := make([]byte, byteSize)
	ReadFull(&conn, data)

	out := new(ring.Poly)
	out.UnmarshalBinary(data)

	netObj.UpdateReceiverLog(pid, len(buf)+len(data))

	return out
}

func (netObj *Network) ReceivePolyMat(pid int) [][]ring.Poly {
	conn := netObj.conns[pid]

	//sizes
	buf := make([]byte, 8)
	ReadFull(&conn, buf)
	byteSize := binary.LittleEndian.Uint64(buf)
	sizesData := make([]byte, byteSize)
	ReadFull(&conn, sizesData)

	//poly values
	ReadFull(&conn, buf)
	polyBytes := binary.LittleEndian.Uint64(buf)
	polyData := make([]byte, polyBytes)
	ReadFull(&conn, polyData)

	netObj.UpdateReceiverLog(pid, 2*len(buf)+len(polyData)+len(sizesData))

	return UnmarshalPolyMat(sizesData, polyData)
}

// func SendPolyMat(poly [][]*ring.Poly, pid int) {
// 	conn := conns[pid]
// 	//marshal poly
// 	sizes, bytes := MarshalPolyMat(poly)

// 	//send over sizes
// 	sbuf := make([]byte, 8)

// 	binary.LittleEndian.PutUint64(sbuf, uint64(len(sizes)))
// 	WriteFull(&conn, sbuf)
// 	WriteFull(&conn, sizes)

// 	binary.LittleEndian.PutUint64(sbuf, uint64(len(bytes)))
// 	WriteFull(&conn, sbuf)
// 	WriteFull(&conn, bytes)

// 	netObj.UpdateSenderLog(pid, len(bytes))

// }

// func ReceivePolyMat(pid int) [][]*ring.Poly {
// 	conn := conns[pid]

// 	//sizes
// 	buf := make([]byte, 8)
// 	ReadFull(&conn, buf)
// 	byteSize := binary.LittleEndian.Uint64(buf)
// 	sizesData := make([]byte, byteSize)
// 	ReadFull(&conn, sizesData)

// 	//poly values
// 	ReadFull(&conn, buf)
// 	polyBytes := binary.LittleEndian.Uint64(buf)
// 	polyData := make([]byte, polyBytes)
// 	ReadFull(&conn, polyData)

// 	netObj.UpdateReceiverLog(pid, len(sizesData))
// 	netObj.UpdateReceiverLog(pid, len(polyData))

// 	return UnmarshalPolyMat(sizesData, polyData)
// }

//SEND AND RECEIVE FIELD VALUES

//ReceiveRMat receives matrix from party p
func (netObj *Network) ReceiveRMat(rtype mpc_core.RElem, n, m, p int) mpc_core.RMat {
	conn := netObj.conns[p]

	buf := make([]byte, 4)
	ReadFull(&conn, buf)

	byteSize := binary.LittleEndian.Uint32(buf)
	data := make([]byte, byteSize)
	ReadFull(&conn, data)

	out := mpc_core.InitRMat(rtype.Zero(), n, m)
	if byteSize != out.NumBytes() {
		panic(fmt.Sprint("Data received:", byteSize, "bytes, Expected", out.NumBytes(), "bytes"))
	}
	out.UnmarshalBinary(data)

	netObj.UpdateReceiverLog(p, len(buf)+len(data))

	return out
}

func (netObj *Network) ReceiveRVec(rtype mpc_core.RElem, n, p int) mpc_core.RVec {
	conn := netObj.conns[p]

	buf := make([]byte, 4)
	ReadFull(&conn, buf)

	byteSize := binary.LittleEndian.Uint32(buf)
	data := make([]byte, byteSize)
	ReadFull(&conn, data)

	out := mpc_core.InitRVec(rtype.Zero(), n)
	if byteSize != out.NumBytes() {
		panic(fmt.Sprint("Data received:", byteSize, "bytes, Expected", out.NumBytes(), "bytes"))
	}
	out.UnmarshalBinary(data)

	netObj.UpdateReceiverLog(p, len(buf)+len(data))

	return out
}

func (netObj *Network) ReceiveRElem(rtype mpc_core.RElem, p int) mpc_core.RElem {
	conn := netObj.conns[p]

	buf := make([]byte, rtype.NumBytes())
	ReadFull(&conn, buf)

	netObj.UpdateReceiverLog(p, len(buf))
	return rtype.FromBytes(buf)
}

func (netObj *Network) ExchangeCiphertext(cryptoParams *crypto.CryptoParams, ct *ckks.Ciphertext, all bool) []*ckks.Ciphertext {
	pid := netObj.pid
	aggData := make([]*ckks.Ciphertext, netObj.NumParties)
	for i := 1; i < netObj.NumParties; i++ {
		if !all && pid == 0 {
			return nil
		}
		if i == pid {
			aggData[i] = ct
		} else if i > pid {
			//receive first and then send
			aggData[i] = netObj.ReceiveCiphertext(cryptoParams, i)
			netObj.SendCiphertext(ct, i)
		} else {
			//send and then receive
			netObj.SendCiphertext(ct, i)
			aggData[i] = netObj.ReceiveCiphertext(cryptoParams, i)
		}
	}
	return aggData

}

//ExchangeCipherVector takes in number of cts for vec
func (netObj *Network) ExchangeCipherVector(cryptoParams *crypto.CryptoParams, cv crypto.CipherVector, nct int) []crypto.CipherVector {
	pid := netObj.pid
	aggData := make([]crypto.CipherVector, netObj.NumParties)

	for i := 1; i < netObj.NumParties; i++ {
		if i == pid {
			aggData[i] = cv
		} else if i > pid {
			//receive first and then send
			aggData[i] = netObj.ReceiveCipherVector(cryptoParams, nct, i)
			netObj.SendCipherVector(cv, i)
		} else {
			//send and then receive
			netObj.SendCipherVector(cv, i)
			aggData[i] = netObj.ReceiveCipherVector(cryptoParams, nct, i)
		}
	}
	return aggData
}

//ExchangeCipherMatrix takes in number of vectors = len(cm) and number of cts per vec = len(cm[0])
func (netObj *Network) ExchangeCipherMatrix(cryptoParams *crypto.CryptoParams, cm crypto.CipherMatrix, ncv, nct int) []crypto.CipherMatrix {
	pid := netObj.pid
	aggData := make([]crypto.CipherMatrix, netObj.NumParties)
	for i := 1; i < netObj.NumParties; i++ {
		if i == pid {
			aggData[i] = cm
		} else if i > pid {
			//receive first and then send
			aggData[i] = netObj.ReceiveCipherMatrix(cryptoParams, ncv, nct, i)
			netObj.SendCipherMatrix(cm, i)
		} else {
			//send and then receive
			netObj.SendCipherMatrix(cm, i)
			aggData[i] = netObj.ReceiveCipherMatrix(cryptoParams, ncv, nct, i)
		}
	}
	return aggData
}
