#include "sockets.h"

/*
**********************************************************
函数名：打开网口
参数：sock    网口名
      IP      网口IP
	  Port    网口端口
返回值：true—打开成功，false—打开失败
函数功能：打开对应IP和端口的网口
**********************************************************
*/
bool OpenSocket(SOCKET& sock, const char IP[], const unsigned short Port)
{
	WSADATA wsaData;
	SOCKADDR_IN addrSrv;

	if (!WSAStartup(MAKEWORD(1, 1), &wsaData))
	{
		if ((sock = socket(AF_INET, SOCK_STREAM, 0)) != INVALID_SOCKET)
		{
			addrSrv.sin_addr.S_un.S_addr = inet_addr(IP);
			addrSrv.sin_family = AF_INET;
			addrSrv.sin_port = htons(Port);
			connect(sock, (SOCKADDR*)&addrSrv, sizeof(SOCKADDR));
			return true;
		}
	}
	return false;
}

/*
**********************************************************
函数名：关闭网口
参数：sock    网口名
函数功能：关闭网口名对应的网口
**********************************************************
*/
void CloseSocket(SOCKET& sock)
{
	closesocket(sock);
	WSACleanup();
}