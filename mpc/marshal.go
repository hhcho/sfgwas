package mpc

import (
	"encoding/binary"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"github.com/ldsec/lattigo/v2/ring"
	"go.dedis.ch/onet/v3/log"
)

func MarshalRData(val interface{}) []byte {
	switch t := val.(type) {
	case mpc_core.RElem:
		buf := make([]byte, t.NumBytes())
		t.ToBytes(buf)
		return buf
	case mpc_core.RVec:
		buf, _ := t.MarshalBinary()
		return buf
	case mpc_core.RMat:
		buf, _ := t.MarshalBinary()
		return buf
	default:
		log.LLvl1(time.Now().Format(time.RFC3339), "Cannot marshal unknown type, ", t)
	}
	return nil
}

//Wrappers for marshaling functions defined in crypto

// MarshalCV returns bytes for ct sizes and bytes for each chiphertext: 8 bytes used //TODO: check if 4 bytes are enough
func MarshalCV(cv crypto.CipherVector) ([]byte, []byte) {
	cvBytes, ctSizes, err := cv.MarshalBinary()

	if err != nil {
		panic(err)
	}
	intsize := uint64(8)
	ncts := len(cv)

	sizesbuf := make([]byte, intsize*uint64(ncts))

	offset := uint64(0)
	for i := range ctSizes {
		binary.LittleEndian.PutUint64(sizesbuf[offset:offset+intsize], uint64(ctSizes[i]))
		offset += intsize
	}

	//add len of ciphervector (ie # of ciphertexts)
	return sizesbuf, cvBytes

}

// UnmarshalCipherVector updates cv to have the correct ciphertexts (cv must have the corrext size)
func UnmarshalCV(cryptoParams *crypto.CryptoParams, ncts int, sbytes, ctbytes []byte) crypto.CipherVector {

	intsize := uint64(8)
	offset := uint64(0)
	sizes := make([]int, ncts)
	for i := range sizes {
		sizes[i] = int(binary.LittleEndian.Uint64(sbytes[offset : offset+intsize]))
		offset += intsize
	}

	cv := make(crypto.CipherVector, 1)
	err := (&cv).UnmarshalBinary(cryptoParams, ctbytes, sizes)
	if err != nil {
		panic(err)

	}

	return cv
}

//TODO: dereferenced version --> need to check

func MarshalPolyMat(polys [][]ring.Poly) ([]byte, []byte) {
	sizes := make([]byte, 0)
	polyBytes := make([]byte, 0)

	//add sizes of polys to marshal
	buf := make([]byte, 8)
	binary.LittleEndian.PutUint64(buf, uint64(len(polys)))
	sizes = append(sizes, buf...)
	binary.LittleEndian.PutUint64(buf, uint64(len(polys[0])))
	sizes = append(sizes, buf...)

	polySize := make([]byte, 4)
	for i := range polys {
		for j := range polys[i] {
			poly, _ := polys[i][j].MarshalBinary()

			//stores bytes for polynomial
			polyBytes = append(polyBytes, poly...)

			//number of bytes used for this polynomial
			binary.LittleEndian.PutUint32(polySize, uint32(len(poly)))
			sizes = append(sizes, polySize...)
		}
	}

	return sizes, polyBytes
}

func UnmarshalPolyMat(sizes, bytes []byte) [][]ring.Poly {
	nrows := binary.LittleEndian.Uint64(sizes[0:8])
	ncols := binary.LittleEndian.Uint64(sizes[8:16])

	sizes = sizes[16:]

	polys := make([][]ring.Poly, nrows)

	poffset := uint64(0)
	soffset := uint64(0)
	for i := range polys {
		polys[i] = make([]ring.Poly, ncols)
		for j := range polys[i] {
			polySize := uint64(binary.LittleEndian.Uint32(sizes[soffset : soffset+4]))
			soffset += 4
			polys[i][j] = *new(ring.Poly)
			polys[i][j].UnmarshalBinary(bytes[poffset : poffset+polySize])
			poffset += polySize
		}

	}
	return polys

}

// func MarshalPolyMat(polys [][]*ring.Poly) ([]byte, []byte) {
// 	sizes := make([]byte, 0)
// 	polyBytes := make([]byte, 0)

// 	//add sizes of polys to marshal
// 	buf := make([]byte, 8)
// 	binary.LittleEndian.PutUint64(buf, uint64(len(polys)))
// 	sizes = append(sizes, buf...)
// 	binary.LittleEndian.PutUint64(buf, uint64(len(polys[0])))
// 	sizes = append(sizes, buf...)

// 	polySize := make([]byte, 4)
// 	for i := range polys {
// 		for j := range polys[i] {
// 			poly, _ := polys[i][j].MarshalBinary()

// 			//stores bytes for polynomial
// 			polyBytes = append(polyBytes, poly...)

// 			//number of bytes used for this polynomial
// 			binary.LittleEndian.PutUint32(polySize, uint32(len(poly)))
// 			sizes = append(sizes, polySize...)
// 		}
// 	}

// 	return sizes, polyBytes
// }

// func UnmarshalPolyMat(sizes, bytes []byte) [][]*ring.Poly {
// 	nrows := binary.LittleEndian.Uint64(sizes[0:8])
// 	ncols := binary.LittleEndian.Uint64(sizes[8:16])

// 	sizes = sizes[16:]

// 	polys := make([][]*ring.Poly, nrows)

// 	poffset := uint32(0)
// 	soffset := uint32(0)
// 	for i := range polys {
// 		polys[i] = make([]*ring.Poly, ncols)
// 		for j := range polys[i] {
// 			polySize := binary.LittleEndian.Uint32(sizes[soffset : soffset+4])
// 			soffset += 4
// 			polys[i][j] = new(ring.Poly)
// 			polys[i][j].UnmarshalBinary(bytes[poffset : poffset+polySize])
// 			poffset += polySize
// 		}

// 	}
// 	return polys

// }
