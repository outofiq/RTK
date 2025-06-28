#include<stdio.h>
#include<windows.h>

#pragma comment(lib,"WS2_32.lib")
#pragma warning(disable:4996)


//打开网口
bool OpenSocket(SOCKET& sock, const char IP[], const unsigned short Port);
//关闭网口
void CloseSocket(SOCKET & sock);
