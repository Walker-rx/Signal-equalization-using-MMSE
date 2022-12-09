#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <WinSock2.h>
#include <Ws2tcpip.h>
#include <Windows.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#pragma comment(lib,"ws2_32.lib")  
#include <tchar.h>
#include <string>
#include <mex.h>
using namespace std;

const int InputBufferSize = 8192;
const int remotePort = 8081;
const int commandPort = 1234;
const int PcPort = 8080;     //must be 8080
//const int stopAddr_0 = 0x00810b00;
const int offset = 0x2AF;
//int stopAddr = stopAddr_0 + offset;
const int Frame_len = 1408;
const int minSize = Frame_len / 2;
//const int n_frame = (int)((double)stopAddr  / 8 / Frame_len*64);
const int timeout = 8000;  //ms
int nRcvBuf = 16 * 1024 * 1024;
int startSendReceive(){
	WORD socketVersion = MAKEWORD(2, 2);
	WSADATA wsaData;  // initialize Winsock
	if (WSAStartup(socketVersion, &wsaData) != 0)
	{
		return 0;
	}
	SOCKET serSocket = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
	if (serSocket == INVALID_SOCKET) {
		printf("socket error!");
		return -1;
	}
	sockaddr_in senAddr,command_addr,RecAddr,remoteAddr;
    //控制�??
	command_addr.sin_family = AF_INET;
	command_addr.sin_port = htons(commandPort);
	// command_addr.sin_addr.S_un.S_addr = inet_addr("192.168.0.2");
	InetPton(AF_INET, _T("192.168.0.2"), &command_addr.sin_addr.s_addr);
    //数据�??
    senAddr.sin_family = AF_INET;
	senAddr.sin_port = htons(remotePort);
	//senAddr.sin_addr.S_un.S_addr = inet_addr("192.168.0.2");
	InetPton(AF_INET, _T("192.168.0.2"), &senAddr.sin_addr.s_addr);
    //PC接收�??
    RecAddr.sin_family = AF_INET;
	RecAddr.sin_port = htons(PcPort);
	// RecAddr.sin_addr.S_un.S_addr = inet_addr("192.168.0.3");
	InetPton(AF_INET, _T("192.168.0.3"), &RecAddr.sin_addr.s_addr);
    //receive
    int nAddrLen = sizeof(remoteAddr);

    if ( bind (serSocket, (sockaddr*)&RecAddr, sizeof(RecAddr)) == SOCKET_ERROR) {
		printf("Bind Error!\n");
		closesocket(serSocket);
		return -1;
	}
 	setsockopt(serSocket, SOL_SOCKET, SO_RCVBUF, (char*)&nRcvBuf, sizeof(nRcvBuf));
	int setret = setsockopt(serSocket, SOL_SOCKET, SO_RCVTIMEO, (char*)&timeout, sizeof(timeout));
	//cout << setret << endl;
	int addrLen = sizeof(senAddr);   
	string catalog(".\\data_bin\\");  //***********************catalog catalog catalog
	string txfilename("transmit_test.bin");
    string path_tx = catalog+txfilename;         //*********************address of file

    //发�?�配�??
    short *DataBuf = new short[minSize];
	char *sendData = new char[2 * minSize];
	FILE *pFile_tx;
	pFile_tx = fopen(path_tx.data(), "rb");
	int dataLength = 0;
	int ret = 1;
	int i = 0, j = 0;
	int Fcnt = 0;
	int temp_i = 0;
    // send ddr stop address
	string ddraddrfilename("ddr_stop_addr_test.bin");
	string ddraddrpath = catalog+ddraddrfilename;
	FILE *pFile_addr;
	pFile_addr = fopen(ddraddrpath.data(), "rb");
	int ddraddr_tmp[1];
	fread(ddraddr_tmp, sizeof(int), 1, pFile_addr);
	fclose(pFile_addr);
	char command_data[4];  //ddr stop address
	for(i=0;i<4;i++){
		command_data[i] = (char)(ddraddr_tmp[0]>>(3-i)*8)&0x000000ff;
	}
	int stopAddr_0 = ddraddr_tmp[0];
	int stopAddr = stopAddr_0 + offset;
	int n_frame = (int)((double)stopAddr  / 8 / Frame_len*64)-1;
	sendto(serSocket,command_data,4,0,(sockaddr*)&command_addr,addrLen);
    //send data
    do {
		temp_i++;
		ret = fread(DataBuf, sizeof(short), minSize, pFile_tx);
		dataLength += ret;
		for (i = 0; i < ret; i++) {
			sendData[2 * i+ 1] = (char)(((DataBuf[i]) & (0xff00)) >> 8);
			sendData[2 * i ] = (char)(((DataBuf[i]) & (0x00ff)));
		}
		sendto(serSocket, sendData, Frame_len, 0, (sockaddr*)&senAddr, addrLen);
		Fcnt++;
	} while (ret > 0);
    //cout << "number of tx frames: " <<Fcnt << endl;
	
	delete[] DataBuf;
	delete[] sendData;
	fclose(pFile_tx);

	string rxfilename("receive_test.bin");
	string path_rx = catalog+rxfilename;
	remove(path_rx.data());



	char* recvData = new char[InputBufferSize];
	i = 0;
	double time = 0;
	double counts = 0;
	short *ADCdata = new short[50*1024*1024];
	LARGE_INTEGER nFreq;
	LARGE_INTEGER nBeginTime;
	LARGE_INTEGER nEndTime;
	QueryPerformanceFrequency(&nFreq);
	QueryPerformanceCounter(&nBeginTime);//�??始计�?? 
	FILE *pFile_rx;
	pFile_rx = fopen(path_rx.data(), "ab");			//*********************文件地址
	if (pFile_rx == NULL) {
		printf("can't open!");
		return -1;
	}
	//printf("receiving...\n");
	int tmp_i=0;
	int fcnt = 0;
	for (i = 0; i < n_frame;i++) {
		ret = recvfrom(serSocket, recvData, InputBufferSize, 0, (sockaddr*)&remoteAddr, &nAddrLen);
		if (ret==-1) { 
			printf("timeout!\n");
			break; 
		}
		for (j = 0; j < ret/2; j++) {
			ADCdata[ret/2*tmp_i + j] = ((recvData[2 * j] & 0x003f) << 8) + (recvData[2 * j+1] & 0xff) -((recvData[2*j] & 0x003f) >= 32) * 16384;

		}
		tmp_i++;
		if((i+1) % 10000 == 0){
			fwrite(ADCdata, sizeof(short), 10000*ret/2, pFile_rx);
			tmp_i = 0;
		}

		fcnt++;

	}
	fwrite(ADCdata, sizeof(short), (fcnt % 10000) * 704, pFile_rx);
	fclose(pFile_rx);
	QueryPerformanceCounter(&nEndTime);            //停止计时  
	time = (double)(nEndTime.QuadPart - nBeginTime.QuadPart) / (double)nFreq.QuadPart;//计算程序执行时间单位为s 
	//cout << "recv time: " << time << "s" << endl;
	//printf("Frame length : %d bytes\n", ret);
	//printf("recv frame counts : %d\n", fcnt);
	double rate;
	rate = (double)n_frame / time * 8 * ret;
	//printf("transmission rate is(with saving file) :%f Mbps\n", rate / 1e6);
	//printf("loss ratio : %f%%\n",100*(1-(double)fcnt / n_frame));

	delete[] ADCdata;
	delete[] recvData;
	closesocket(serSocket);
	WSACleanup();
    return 0;
}
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ){
    startSendReceive();
	return;
}