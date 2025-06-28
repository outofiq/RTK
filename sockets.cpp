#include "sockets.h"

/*
**********************************************************
��������������
������sock    ������
      IP      ����IP
	  Port    ���ڶ˿�
����ֵ��true���򿪳ɹ���false����ʧ��
�������ܣ��򿪶�ӦIP�Ͷ˿ڵ�����
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
���������ر�����
������sock    ������
�������ܣ��ر���������Ӧ������
**********************************************************
*/
void CloseSocket(SOCKET& sock)
{
	closesocket(sock);
	WSACleanup();
}