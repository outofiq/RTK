#include "header.h"

using namespace std;

/*
**************************************************************************************
�������������ļ���ȡ����
������filename �����ļ�·��������
	  set      �洢���ò����Ľṹ�����
����ֵ��true----��ȡ�����ļ��ɹ���false---��ȡ�����ļ�ʧ��
�������ܣ���ȡ�����ļ��������Ӧ��������ȷ��������ʽ�����ļ���ʽ����������ʽ����ϵͳ����
**************************************************************************************
*/
bool ReadCFG(const char cfgname[],CFGINFO_t* cfg)
{
	ifstream infile(cfgname);
	if (!infile.is_open())
	{
		cerr << "open error!" << endl;
		return false;
	}

	string header;
	getline(infile, header);

	string line;
	getline(infile, line);
	istringstream iss(line);
	iss >> cfg->datatype
		>> cfg->systype
		>> cfg->BasPos_x
		>> cfg->BasPos_y
		>> cfg->BasPos_z
		>> cfg->BasIP
		>> cfg->BasPort
		>> cfg->RovIP
		>> cfg->RovPort
		>> cfg->BasInfo
		>> cfg->RovInfo
		>> cfg->RTKType
		>> cfg->dx
		>> cfg->dy
		>> cfg->dz
		>> cfg->dE
		>> cfg->dN
		>> cfg->dU;
	infile.close();
	return true;
}