package mpc

import (
	// "bufio"

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
