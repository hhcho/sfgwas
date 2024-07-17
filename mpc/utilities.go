package mpc

import (
	// "bufio"

	"github.com/hhcho/sfgwas/crypto"
	"strconv"
	"strings"
)

// for processing terminal input => to get pid and ip addresses and port
func split(s string) (int, string, string) {
	info := strings.SplitN(s, ":", 3)
	// fmt.Println(info)
	id, _ := strconv.Atoi(info[0])
	return id, info[1], info[2]
}

func checkError(e error) {
	if e != nil {
		panic(e)
	}
}

func Debug(cryptoParams *crypto.CryptoParams, mpcObj *MPC, matrix crypto.CipherMatrix, indexes []int, name string,
	debug bool, numParties int) {

	if debug {
		for party := 1; party <= numParties; party++ {
			if len(indexes) > 1 {
				SaveMatrixToFileWithPrint(cryptoParams, mpcObj, matrix, indexes[party],
					party, name, true)
			} else {
				SaveMatrixToFileWithPrintIndex(cryptoParams, mpcObj, matrix, indexes[0],
					party, name, true, 0, indexes[1])
			}
		}
	}
}
