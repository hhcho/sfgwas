package mpc

import (
	// "bufio"

	"fmt"
	"net"
	"os"
	"strconv"
	"sync"
	"time"

	"github.com/aead/chacha20/chacha"
	"github.com/hhcho/frand"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/dckks"
	"github.com/ldsec/lattigo/v2/ring"
	"go.dedis.ch/onet/v3/log"
)

type Network struct {
	pid        int
	hubPid     int
	NumParties int

	crpGen       *ring.UniformSampler
	dckksContext *dckks.Context
	Rand         *Random

	// Paralellization: socket between every machine, thread pair

	conns     map[int]net.Conn
	listeners map[int]net.Listener

	// Variables to keep track of the bytes that are sent/received

	SentBytes     map[int]uint64
	ReceivedBytes map[int]uint64
	commSent      map[int]int
	commReceived  map[int]int

	loggingActive bool
}

func (netObj *Network) EnableLogging() {
	netObj.loggingActive = true
}

func (netObj *Network) DisableLogging() {
	netObj.loggingActive = false
}

func (netObj *Network) UpdateSenderLog(toPid int, nbytes int) {
	if netObj.loggingActive {
		netObj.SentBytes[toPid] += uint64(nbytes)
		netObj.commSent[toPid]++
	}
}

func (netObj *Network) UpdateReceiverLog(fromPid int, nbytes int) {
	if netObj.loggingActive {
		netObj.ReceivedBytes[fromPid] += uint64(nbytes)
		netObj.commReceived[fromPid]++
	}
}

func (netObj *Network) ResetNetworkLog() {
	for key := range netObj.SentBytes {
		netObj.SentBytes[key] = 0
	}
	for key := range netObj.ReceivedBytes {
		netObj.ReceivedBytes[key] = 0
	}
}

func (netObjs ParallelNetworks) ResetNetworkLog() {
	for i := range netObjs {
		netObjs[i].ResetNetworkLog()
	}
}

func (netObjs ParallelNetworks) PrintNetworkLog() {
	sentBytes := make(map[int]uint64)
	receivedBytes := make(map[int]uint64)

	for i := range netObjs {
		for key, value := range netObjs[i].SentBytes {
			sentBytes[key] += value
		}
		for key, value := range netObjs[i].ReceivedBytes {
			receivedBytes[key] += value
		}
	}

	fmt.Println("Network log for party", netObjs[0].pid)
	for key, value := range sentBytes {
		fmt.Println(value, "bytes to party", key)
	}
	for key, value := range receivedBytes {
		fmt.Println(value, "bytes from party", key)
	}

}

func (netObj *Network) PrintNetworkLog() {
	fmt.Println("Network log for party", netObj.pid)
	for key, value := range netObj.SentBytes {
		fmt.Println(value, "bytes to party", key)
	}
	for key, value := range netObj.ReceivedBytes {
		fmt.Println(value, "bytes from party", key)
	}
}

type Server struct {
	IpAddr string
	Ports  map[string]string
}

func pidString(pid int) string {
	return fmt.Sprintf("party%d", pid)
}

/* Communication set up for parallelization */

// InitNetworkParalelization creates communication channels for parallelization: channel btw each pair of machines for every thread
func InitCommunication(bindingIP string, servers map[string]Server, pid, nparties, numThreads int, sharedKeysPath string) []*Network {
	network := make([]*Network, numThreads)

	if bindingIP == "" { // Set default
		bindingIP = "0.0.0.0"
	}

	// Add sync to make sure communication for current thread is set
	var wg sync.WaitGroup
	for thread := range network {
		wg.Add(1)

		go func(thread int) {
			defer wg.Done()

			// Create network layer for this thread
			network[thread] = initNetworkForThread(bindingIP, servers, pid, nparties, thread)

			fmt.Printf("Network for thread %v complete. ", thread)
		}(thread)
	}

	wg.Wait()

	// update threads initialize using different outputs of the same seed
	InitializeParallelPRG(sharedKeysPath, network, pid, nparties)

	return network

}

func InitializeParallelPRG(sharedKeysPath string, network []*Network, pid int, nparties int) {
	if sharedKeysPath == "" {
		log.LLvl1("Warning: shared_keys_path not set in config. Falling back on deterministic keys (not secure).")
	}
	randMaster := InitializePRG(pid, nparties, sharedKeysPath)
	for i := range network {
		for j := -1; j < nparties; j++ {
			seed := make([]byte, chacha.KeySize)
			randMaster.SwitchPRG(j)
			randMaster.RandRead(seed)
			randMaster.RestorePRG()
			network[i].Rand.prgTable[j] = frand.NewCustom(seed, bufferSize, 20)
		}
		network[i].Rand.curPRG = network[i].Rand.prgTable[pid]
		network[i].Rand.pid = pid
	}
}

func initNetworkForThread(bindingIP string, servers map[string]Server, pid int, nparties, thread int) *Network {
	conns := make(map[int]net.Conn)
	listeners := make(map[int]net.Listener)

	fmt.Println("Initializing network for thread: ", thread)

	for other := 0; other < nparties; other++ {
		if other == pid {
			continue
		}

		if other < pid { // Act as a client
			ip := servers[pidString(other)].IpAddr
			if ip == servers[pidString(pid)].IpAddr { // If self, set to localhost
				ip = "127.0.0.1"
			}

			portInt, error := strconv.Atoi(servers[pidString(other)].Ports[pidString(pid)])
			checkError(error)

			// Need different port for each pair of thread between each pair of parties
			// Assumes that port info in config is given such that port + thread ID does not collide
			port := strconv.Itoa(portInt + thread)

			conns[other] = Connect(ip, port)
			fmt.Println("Connected to socket, listening to party", other)

		} else { // Act as a server
			ip := bindingIP
			portInt, error := strconv.Atoi(servers[pidString(pid)].Ports[pidString(other)])
			checkError(error)

			port := strconv.Itoa(portInt + thread)
			c, l := OpenChannel(ip, port)
			if c == nil || l == nil {
				panic("Opening communication channels was unsuccessful")
			}
			conns[other] = c
			listeners[other] = l

			fmt.Println("Check port: opened and listening to party", other)
		}

	}

	return &Network{
		conns:         conns,
		listeners:     listeners,
		NumParties:    nparties,
		pid:           pid,
		hubPid:        1,
		Rand:          nil,
		SentBytes:     make(map[int]uint64),
		ReceivedBytes: make(map[int]uint64),
		commSent:      make(map[int]int),
		commReceived:  make(map[int]int),
		loggingActive: true,
	}

}

func (netObj *Network) SetMHEParams(params *ckks.Parameters) {
	dckksContext := dckks.NewContext(params)

	seed := make([]byte, chacha.KeySize)
	netObj.Rand.SwitchPRG(-1)
	netObj.Rand.RandRead(seed)
	netObj.Rand.RestorePRG()

	seedPrng := frand.NewCustom(seed, bufferSize, 20)

	crpGen := ring.NewUniformSamplerWithBasePrng(seedPrng, dckksContext.RingQP)

	netObj.crpGen = crpGen
	netObj.dckksContext = dckksContext
}

/* Reading and writing from channel */

func WriteFull(conn *net.Conn, buf []byte) {
	shift := 0
	remaining := len(buf)
	for {
		sent, err := (*conn).Write(buf[shift:])
		if sent == remaining {
			break
		} else if sent < remaining {
			shift += sent
			remaining -= sent
		}

		if err != nil {
			panic(err)
		}
	}
}

func ReadFull(conn *net.Conn, buf []byte) {
	shift := 0
	remaining := len(buf)
	for {
		received, err := (*conn).Read(buf[shift:])
		if received == remaining {
			break
		} else if received < remaining {
			shift += received
			remaining -= received
		}

		if err != nil {
			panic(err)
		}
	}
}

func ReadFullWithErr(conn *net.Conn, buf []byte) error {
	shift := 0
	remaining := len(buf)
	for {
		received, err := (*conn).Read(buf[shift:])
		if received == remaining {
			break
		} else if received < remaining {
			shift += received
			remaining -= received
		}

		if err != nil {
			return err
		}
	}
	return nil
}

//OpenChannel opens channel at specificed ip address and port and returns channel (server side, connection for client to listen to)
func OpenChannel(ip, port string) (net.Conn, net.Listener) {
	l, err := establishConn(ip, port)
	checkError(err)
	conn, err := listen(l)
	checkError(err)
	fmt.Println("Successfully opened channel at port " + port)

	//return the channel specific to party pair
	return conn, l
}

func establishConn(ip, port string) (net.Listener, error) {
	addr := ":" + port
	if ip != "" {
		addr = ip + addr

		fmt.Println("Opening socket on: ", addr)
	}

	l, err := net.Listen("tcp", addr)
	if err != nil {
		fmt.Println(err)
		return nil, err
	}

	return l, nil
}

func listen(l net.Listener) (net.Conn, error) {
	c, err := l.Accept()
	if err != nil {
		// fmt.Println(err)
		return nil, err
	}
	return c, nil
}

//CloseChannel closes connection
func CloseChannel(conn net.Conn) {
	err := conn.Close()
	checkError(err)
}

func (netObj *Network) CloseAll() {
	for _, c := range netObj.conns {
		CloseChannel(c)
	}
	for _, l := range netObj.listeners {
		l.Close()
	}
}

// Connect to "server", given the ip address and port
func Connect(ip, port string) net.Conn {
	addr := ip + ":" + port

	// Re-try connecting to server
	retrySchedule := make([]time.Duration, 100)
	for i := range retrySchedule {
		// retrySchedule[i] = time.Duration((i+1)*5) * time.Second
		retrySchedule[i] = time.Duration(5) * time.Second
	}

	var c net.Conn
	var err error
	for _, retry := range retrySchedule {
		c, err = net.Dial("tcp", addr)

		if err == nil {
			fmt.Println("Successfully connected to " + addr)
			break
		}

		fmt.Printf("Connection failed. Error: %+v\n", err)
		fmt.Printf("Retrying in... %v\n", retry)
		time.Sleep(retry)
	}

	return c
}

func SaveBytesToFile(b []byte, filename string) {
	f, err := os.Create(filename)

	if err != nil {
		panic(err)
	}

	defer f.Close()

	f.Write(b)

	f.Sync()
}

func (netObj *Network) GetConn(to int, threadNum int) net.Conn {
	return netObj.conns[to]
}

func (netObj *Network) SetPid(p int) {
	netObj.pid = p
}

func (netObj *Network) GetPid() int {
	return netObj.pid
}

func (netObj *Network) SetHubPid(p int) {
	netObj.hubPid = p
}

func (netObj *Network) GetHubPid() int {
	return netObj.hubPid
}

func (netObj *Network) SetNParty(np int) {
	netObj.NumParties = np
}

func (netObj *Network) GetNParty() int {
	return netObj.NumParties
}

func (netObj *Network) GetCRPGen() *ring.UniformSampler {
	return netObj.crpGen
}
