package gwas

import (
	"bufio"
	"encoding/binary"
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"strconv"
	"strings"

	"github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"gonum.org/v1/gonum/mat"
)

func fileExists(filename string) bool {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}

// Reads in a binary file containing 6 vectors of length m (# of SNPs):
// ref allele count (AC), alt AC, hom-ref genotype count (GC), het GC,
// hom-alt GC, missing sample count.
// Each value is encoded as uint32 in little endian format
func ReadGenoStatsFromFile(filename string, m int) (ac, gc [][]uint32, miss []uint32) {
	nstats := 6

	file, err := os.Open(filename)
	defer file.Close()
	if err != nil {
		log.Fatal(err)
	}

	reader := bufio.NewReader(file)

	out := make([][]uint32, nstats)
	buf := make([]byte, 4*m) // 4 bytes per number
	for s := 0; s < nstats; s++ {
		out[s] = make([]uint32, m)

		if _, err := io.ReadFull(reader, buf); err != nil {
			log.Fatal(err)
		}

		for i := range out[s] {
			out[s][i] = binary.LittleEndian.Uint32(buf[4*i:])
		}
	}

	gc = out[:3]

	// TODO: remove allele counts from file as they are redundant
	for i := range out[3] {
		out[3][i] = out[1][i] + 2*out[0][i]
		out[4][i] = out[1][i] + 2*out[2][i]
	}
	ac = out[3:5]

	miss = out[5]

	return
}

func readFilterFromFile(filename string, n int, isBinary bool) []bool {
	file, err := os.Open(filename)
	defer file.Close()
	if err != nil {
		log.Fatal(err)
	}

	reader := bufio.NewReader(file)

	filter := make([]bool, n)

	if isBinary {

		buf := make([]byte, len(filter))
		if _, err := io.ReadFull(reader, buf); err != nil {
			log.Fatal(err)
		}

		for i := range filter {
			filter[i] = buf[i] != 0
		}

	} else {

		for i := range filter {
			line, err := reader.ReadString('\n')
			if err != nil {
				log.Fatal(err)
			}

			filter[i] = line[0] != '0'
		}

	}

	return filter
}

func writeFilterToFile(filename string, filter []bool, isBinary bool) {
	file, err := os.Create(filename)
	defer file.Close()
	if err != nil {
		log.Fatal(err)
	}

	writer := bufio.NewWriter(file)

	if isBinary {
		buf := make([]byte, len(filter))
		for i := range filter {
			if filter[i] {
				buf[i] = byte(1)
			} else {
				buf[i] = byte(0)
			}
		}
		writer.Write(buf)
	} else {
		for i := range filter {
			if filter[i] {
				writer.WriteString("1\n")
			} else {
				writer.WriteString("0\n")
			}
		}
	}

	writer.Flush()
}

func FilterMatrixFilePgen(pgenPrefix string, nrows, ncols int, rowFiltFile, colNamesFile string, colStartPos int, colFilt []bool, outputFile string) {
	colFiltFile := outputFile + ".colFilter.bin"

	writeFilterToFile(colFiltFile, colFilt, true)

	cmd := exec.Command("/bin/sh", "scripts/filterMatrixPgen.sh", pgenPrefix, strconv.Itoa(nrows), strconv.Itoa(ncols), rowFiltFile, colFiltFile, colNamesFile, strconv.Itoa(colStartPos), outputFile)
	cout, e := cmd.CombinedOutput()
	fmt.Print(string(cout))
	if e != nil {
		log.Fatalln(e)
	}
}

func FilterMatrixFile(inputFile string, nrows, ncols int, rowFilt, colFilt []bool, outputFile string) {
	rowFiltFile := inputFile + ".pcaRowFilter.bin"
	colFiltFile := inputFile + ".pcaColFilter.bin"

	writeFilterToFile(rowFiltFile, rowFilt, true)
	writeFilterToFile(colFiltFile, colFilt, true)

	cmd := exec.Command("/bin/sh", "scripts/filterMatrix.sh", inputFile, strconv.Itoa(nrows), strconv.Itoa(ncols), rowFiltFile, colFiltFile, outputFile)
	cout, e := cmd.CombinedOutput()
	fmt.Print(string(cout))
	if e != nil {
		log.Fatalln(e)
	}
}

func TransposeMatrixFile(inputFile string, nrows, ncols int, outputFile string) {
	cmd := exec.Command("/bin/sh", "scripts/transposeMatrix.sh", inputFile, strconv.Itoa(nrows), strconv.Itoa(ncols), outputFile)
	cout, e := cmd.CombinedOutput()
	fmt.Print(string(cout))
	if e != nil {
		log.Fatalln(e)
	}
}

func MergeBlockFiles(inputBlockFilePrefix string, nrows int, ncolsPerBlock []int, outputFile string) {
	datFile := inputBlockFilePrefix + ".merge.sizes.txt"

	// Write dimensions to file
	file, err := os.Create(datFile)
	defer file.Close()
	if err != nil {
		log.Fatal(err)
	}

	writer := bufio.NewWriter(file)
	for _, v := range ncolsPerBlock {
		writer.WriteString(strconv.Itoa(v) + "\n")
	}
	writer.Flush()

	cmd := exec.Command("/bin/sh", "scripts/mergeMatrices.sh", inputBlockFilePrefix, strconv.Itoa(nrows), datFile, outputFile)
	cout, e := cmd.CombinedOutput()
	fmt.Print(string(cout))
	if e != nil {
		log.Fatalln(e)
	}
}

func matPrint(X mat.Matrix) {
	fa := mat.Formatted(X, mat.Prefix(""), mat.Squeeze())
	fmt.Printf("%v\n", fa)
}
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
func Sum(a []int) int {
	res := 0
	for i := range a {
		res += a[i]
	}
	return res
}
func SumF(a []float64) float64 {
	res := 0.0
	for i := range a {
		res += a[i]
	}
	return res
}

func ScaleF(a []float64, scale float64) []float64 {
	res := make([]float64, len(a))
	for i := range a {
		res[i] = a[i] * scale
	}
	return res
}

func SumBool(a []bool) int {
	res := 0
	for i := range a {
		if a[i] {
			res++
		}
	}
	return res
}

func Ones(n int) []int {
	a := make([]int, n)
	for i := range a {
		a[i] = 1
	}
	return a
}

func OnesBool(n int) []bool {
	a := make([]bool, n)
	for i := range a {
		a[i] = true
	}
	return a
}

func LoadCacheFromFile(cps *crypto.CryptoParams, filename string) crypto.CipherMatrix {
	delim := ','

	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}

	defer f.Close()

	// data := make([]float64, 0)

	// See how we can read strings and parse them
	c := csv.NewReader(f)
	c.Comma = delim
	text, err := c.ReadAll()

	columns := c.FieldsPerRecord
	lines := len(text)

	data := make([][]float64, lines)

	for i := 0; i < lines; i++ {
		data[i] = make([]float64, columns)

		for j := 0; j < columns; j++ {
			data[i][j], err = strconv.ParseFloat(text[i][j], 64)
		}
	}

	res, _, _, _ := crypto.EncryptFloatMatrixRow(cps, data)
	return res
}

func LoadMatrixFromFileFloat(filename string, delim rune) [][]float64 {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}

	defer f.Close()

	c := csv.NewReader(f)
	c.Comma = delim
	text, err := c.ReadAll()

	columns := c.FieldsPerRecord
	lines := len(text)

	data := make([][]float64, lines)
	for i := 0; i < lines; i++ {
		data[i] = make([]float64, columns)
		for j := 0; j < columns; j++ {
			data[i][j], err = strconv.ParseFloat(text[i][j], 64)
		}
	}

	return data
}

func LoadMatrixFromFile(filename string, delim rune) *mat.Dense {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}

	defer f.Close()

	c := csv.NewReader(f)
	c.Comma = delim
	text, err := c.ReadAll()

	columns := c.FieldsPerRecord
	lines := len(text)

	data := make([]float64, columns*lines)

	for i := 0; i < lines; i++ {
		for j := 0; j < columns; j++ {
			data[i*columns+j], err = strconv.ParseFloat(text[i][j], 64)
		}
	}

	return mat.NewDense(lines, columns, data)
}

func LoadSNPPositionFile(filename string, delim rune) []uint64 {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}

	defer f.Close()

	c := csv.NewReader(f)
	c.Comma = delim
	text, err := c.ReadAll()

	lines := len(text)

	data := make([]uint64, lines)

	for i := 0; i < lines; i++ {
		chrom, err := strconv.ParseUint(text[i][0], 10, 64)
		if err != nil {
			panic(err)
		}

		pos, err := strconv.ParseUint(text[i][1], 10, 64)
		if err != nil {
			panic(err)
		}

		data[i] = chrom*1e9 + pos
	}

	return data
}

func SaveMatrixToFile(cps *crypto.CryptoParams, mpcObj *mpc.MPC, cm crypto.CipherMatrix, nElemCol int, sourcePid int, filename string) {
	pid := mpcObj.GetPid()
	if pid == 0 {
		return
	}

	pm, _ := mpcObj.Network.CollectiveDecryptMat(cps, cm, sourcePid)

	M := mat.NewDense(len(cm), nElemCol, nil)
	for i := range pm {
		M.SetRow(i, crypto.DecodeFloatVector(cps, pm[i])[:nElemCol])
	}

	if pid == sourcePid || sourcePid < 0 {

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

func SaveFloatMatrixToFile(filename string, x [][]float64) {
	file, err := os.Create(filename)
	defer file.Close()
	if err != nil {
		log.Fatal(err)
	}

	writer := bufio.NewWriter(file)

	for j := range x[0] {
		for i := range x {
			writer.WriteString(fmt.Sprintf("%.6e", x[i][j]))
			if i < len(x)-1 {
				writer.WriteString(",")
			}
		}
		writer.WriteString("\n")
	}

	writer.Flush()
}

func SaveFloatVectorToFile(filename string, x []float64) {
	file, err := os.Create(filename)
	defer file.Close()
	if err != nil {
		log.Fatal(err)
	}

	writer := bufio.NewWriter(file)

	for i := range x {
		writer.WriteString(fmt.Sprintf("%.6e\n", x[i]))
	}

	writer.Flush()
}

func LoadFloatVectorFromFile(filename string, n int) []float64 {
	file, err := os.Open(filename)
	defer file.Close()
	if err != nil {
		log.Fatal(err)
	}

	reader := bufio.NewReader(file)

	out := make([]float64, n)

	for i := range out {
		line, err := reader.ReadString('\n')
		if err != nil {
			log.Fatal(err)
		}

		_, err = fmt.Sscanf(line, "%f", &out[i])
		if err != nil {
			log.Fatal(err)
		}
	}

	return out
}

func SaveIntVectorToFile(filename string, x []int) {
	file, err := os.Create(filename)
	defer file.Close()
	if err != nil {
		log.Fatal(err)
	}

	writer := bufio.NewWriter(file)

	for i := range x {
		writer.WriteString(fmt.Sprintf("%d\n", x[i]))
	}

	writer.Flush()
}
