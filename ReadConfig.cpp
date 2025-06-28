#include "header.h"

using namespace std;

/*
**************************************************************************************
函数名：配置文件读取函数
参数：filename 配置文件路径及名称
	  set      存储配置参数的结构体变量
返回值：true----读取配置文件成功，false---读取配置文件失败
函数功能：读取配置文件，获得相应参数，以确定数据形式（是文件形式还是网口形式）、系统类别等
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