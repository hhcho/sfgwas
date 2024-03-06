package crypto

import (
	"bufio"
	"encoding/binary"
	"io"
	"os"
	"unsafe"

	"go.dedis.ch/onet/v3/log"
)

func Max(x, y int) int {
	if x <= y {
		return y
	}
	return x
}

func Min(a int, b int) int {
	if a > b {
		return b
	}
	return a
}
func Mod(n int, modulus int) int {
	n = n % modulus
	if n < 0 {
		n = n + modulus
	}
	return n
}

// MarshalCiphermatrix returns byte array corresponding to ciphertext sizes (int array) and byte array corresponding to marshaling
func MarshalCM(cm CipherMatrix) ([]byte, []byte) {
	cmBytes, ctSizes, err := cm.MarshalBinary()
	if err != nil {
		panic(err)
	}

	r, c := len(cm), len(cm[0])
	intsize := uint64(unsafe.Sizeof(ctSizes[0][0]))
	sizesbuf := make([]byte, intsize*uint64(r)*uint64(c))

	offset := uint64(0)
	for i := range ctSizes {
		for j := range ctSizes[i] {
			binary.LittleEndian.PutUint64(sizesbuf[offset:offset+intsize], uint64(ctSizes[i][j]))
			offset += intsize
		}

	}

	return sizesbuf, cmBytes

}

// MarshalCiphermatrix returns byte array corresponding to ciphertext sizes (int array) and byte array corresponding to marshaling
func UnmarshalCM(cryptoParams *CryptoParams, r, c int, sbytes, ctbytes []byte) CipherMatrix {
	intsize := uint64(8)
	offset := uint64(0)
	sizes := make([][]int, r)
	for i := range sizes {
		sizes[i] = make([]int, c)
		for j := range sizes[i] {
			sizes[i][j] = int(binary.LittleEndian.Uint64(sbytes[offset:]))
			offset += intsize

		}
	}

	cm := make(CipherMatrix, 1)
	err := (&cm).UnmarshalBinary(cryptoParams, ctbytes, sizes)
	if err != nil {
		panic(err)

	}

	return cm
}

func SaveCipherMatrixToFile(cps *CryptoParams, cm CipherMatrix, filename string) {
	file, err := os.Create(filename)
	defer file.Close()
	if err != nil {
		log.Fatal(err)
	}

	writer := bufio.NewWriter(file)

	sbytes, cmbytes := MarshalCM(cm)

	nrbuf := make([]byte, 4)
	ncbuf := make([]byte, 4)
	binary.LittleEndian.PutUint32(nrbuf, uint32(len(cm)))
	binary.LittleEndian.PutUint32(ncbuf, uint32(len(cm[0])))

	sbuf := make([]byte, 8) //TODO: see if 4 bytes is enough
	cmbuf := make([]byte, 8)
	binary.LittleEndian.PutUint64(sbuf, uint64(len(sbytes)))
	binary.LittleEndian.PutUint64(cmbuf, uint64(len(cmbytes)))

	writer.Write(nrbuf)
	writer.Write(ncbuf)
	writer.Write(sbuf)
	writer.Write(sbytes)
	writer.Write(cmbuf)
	writer.Write(cmbytes)

	writer.Flush()
}

func LoadCipherMatrixFromFile(cps *CryptoParams, filename string) CipherMatrix {
	file, err := os.Open(filename)
	defer file.Close()
	if err != nil {
		log.Fatal(err)
	}

	reader := bufio.NewReader(file)

	ibuf := make([]byte, 4)
	io.ReadFull(reader, ibuf)
	nrows := int(binary.LittleEndian.Uint32(ibuf))
	io.ReadFull(reader, ibuf)
	numCtxPerRow := int(binary.LittleEndian.Uint32(ibuf))

	sbuf := make([]byte, 8)
	io.ReadFull(reader, sbuf)
	sbyteSize := binary.LittleEndian.Uint64(sbuf)
	sdata := make([]byte, sbyteSize)
	io.ReadFull(reader, sdata)

	cmbuf := make([]byte, 8)
	io.ReadFull(reader, cmbuf)
	cbyteSize := binary.LittleEndian.Uint64(cmbuf)
	cdata := make([]byte, cbyteSize)
	io.ReadFull(reader, cdata)

	return UnmarshalCM(cps, nrows, numCtxPerRow, sdata, cdata)
}

//
// // *** FILE MANAGEMENT
//
// // SaveModel saves weights/model to a binary file
// func SaveModel(file string, model CipherMatrix) error {
// 	binaryData, sizes, err := model.MarshalBinary()
//
// 	if err != nil {
// 		return err
// 	}
//
// 	data := make([]byte, binary.MaxVarintLen64*(len(sizes)*(1+len(sizes[0]))+1))
//
// 	index := 0
// 	// Write the int array
// 	index += binary.PutUvarint(data, uint64(len(sizes)))
// 	for _, elt := range sizes {
// 		index += binary.PutUvarint(data[index:], uint64(len(elt)))
// 		for _, e := range elt {
// 			index += binary.PutUvarint(data[index:], uint64(e))
// 		}
// 	}
//
// 	data = data[:index]
// 	data = append(data, binaryData...)
//
// 	return ioutil.WriteFile(file, data, os.ModePerm)
// }
//
// //
// // // LoadModel loads weights/model to a binary file
// func LoadModel(params *CryptoParams, file string) (CipherMatrix, error) {
// 	bts, err := ioutil.ReadFile(file)
// 	if err != nil {
// 		return nil, fmt.Errorf("reading file: %w", err)
// 	}
// 	reader := bytes.NewReader(bts)
//
// 	mainSize, err := binary.ReadUvarint(reader)
// 	if err != nil {
// 		return nil, fmt.Errorf("reading uint: %w", err)
// 	}
//
// 	sizes := make([][]int, mainSize)
// 	for i := range sizes {
// 		childSize, err := binary.ReadUvarint(reader)
// 		if err != nil {
// 			return nil, fmt.Errorf("reading uint: %w", err)
// 		}
//
// 		sizes[i] = make([]int, childSize)
//
// 		for j := range sizes[i] {
// 			size, err := binary.ReadUvarint(reader)
// 			if err != nil {
// 				return nil, fmt.Errorf("reading uint: %w", err)
// 			}
//
// 			sizes[i][j] = int(size)
// 		}
// 	}
//
// 	var mat = make(CipherMatrix, 0)
//
// 	f := make([]byte, reader.Len())
// 	_, err = reader.Read(f)
// 	if err != nil {
// 		return nil, fmt.Errorf("reading file block: %w", err)
// 	}
//
// 	err = mat.UnmarshalBinary(params, f, sizes)
// 	if err != nil {
// 		return nil, fmt.Errorf("unmarshalling weights: %w", err)
// 	}
//
// 	return mat, nil
// }
