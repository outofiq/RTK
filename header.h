#include "timeTrans.h"
#include "coordinateTrans.h"
#include <wchar.h>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include "Matrix.h"
#include "sockets.h"

extern "C"
{
#include "LAMBDA.h";
}

using namespace std;
using namespace Matrix_hx;

#pragma once
/* OEM basic info */
#define OEM7SYNC1       0xAA	     /* oem7/6/4 message start sync code 1 */
#define OEM7SYNC2       0x44		 /* oem7/6/4 message start sync code 2 */
#define OEM7SYNC3       0x12		 /* oem7/6/4 message start sync code 3 */
#define OEM7HLEN        28			 /* oem7/6/4 message header length (bytes) */
#define POLYCRC32       0xEDB88320u

/* message IDs */
#define ID_RANGE		43           /* oem7/6/4 range measurement */
#define ID_BDSEPHEM		1696         /* oem7/6 decoded bds ephemeris */
#define ID_GPSEPHEM		7		     /* oem7 decoded GPS L1 C/A ephemerides */
#define ID_BESTPOS		42		     /* oem7 decoded best position */

#define GPST_BDT		14.0	     /* Difference between GPS time and Beidou time[s] */
#define MAXCHANNUM      36
#define MAXGPSPRN       32
#define MAXGEOPRN       5            /* 最大的GEO卫星号 */ //不完善，BDS-3的GEO被排除在外了
#define MAXBDSPRN       63
#define MAXRAWLEN       40960
#define MAXBUFF	        20480

/* math common uints */
#define PI          3.1415926535898   /* PI */  
#define D2R			(PI/180.0)        /* radians to degree */
#define R2D         (180.0/PI)        /* degree to radians */
#define CLIGHT		299792458.0       /* Speed of light  [m/s]; IAU 1976  */

// GPS卫星信号的一些常量
#define  FG1_GPS  1575.42E6           // L1信号频率 
#define  FG2_GPS  1227.60E6           // L2信号频率
#define  WL1_GPS  (CLIGHT/FG1_GPS)    // L1信号波长
#define  WL2_GPS  (CLIGHT/FG2_GPS)    // L2信号波长

// BDS卫星信号的一些常量
#define  FG1_BDS  1561.098E6          // B1信号的频率 
#define  FG3_BDS  1268.520E6          // B3信号的频率 
#define  WL1_BDS  (CLIGHT/FG1_BDS)    // B1信号波长
#define  WL3_BDS  (CLIGHT/FG3_BDS)    // B3信号波长

#define SECPERHOUR 3600.0             /* Seconds per hour */
#define SECPERDAY  86400.0            /* Seconds per day */
#define SECPERWEEK 604800.0           /* Seconds per week */

/* 卫星系统定义 ---------------------------------------------------------------*/
enum class sys_t { UNKS = 0, GPS, BDS };

//GPS+BDS广播星历
struct gpseph_t 
{
	unsigned short PRN;
	sys_t          Sys;
	GPSTIME        TOC, TOE; //s
	short	       health; //0-健康；1-不健康
	double	       ClkBias, ClkDrift, ClkDriftRate; //钟差 钟速 钟飘
	double	       IODE, IODC; //星历数据龄期 卫星钟数据龄期
	double         TGD1, TGD2; //广播星历播发的时间群延迟
	double	       SqrtA, ecc, M0, OMEGA0, i0, omega, OMEGADot, iDot, DeltaN;
	double	       Crs, Cuc, Cus, Cic, Cis, Crc; //摄动项
	double	       URA;

	gpseph_t()
	{
		PRN = health = 1;
		Sys = sys_t::UNKS;
		ClkBias = ClkDrift = ClkDriftRate = IODE = IODC = TGD1 = TGD2 = 0.0;
		SqrtA = ecc = M0 = OMEGA0 = i0 = omega = OMEGADot = iDot = DeltaN = 0.0;
		Crs = Cuc = Cus = Cic = Cis = Crc = URA = 0.0;
	}
};

//每颗卫星观测数据
struct satobs_t 
{
	sys_t    Sys;
	short    Prn;
	//GPS：L1-0，L2-1   BDS：B1I-0，B3I-1
	double   p[2],l[2],d[2]; //伪距 相位 多普勒
	double   CN0[2];         //载噪比 
	bool     Valid;          // 整体（即双频都在内）有效性


	satobs_t()
	{
		Prn = 0;
		Sys = sys_t::UNKS;
		p[0] = p[1] = l[0] = l[1] = d[0] = d[1] = 0.0;
		CN0[0] = CN0[1] = 0.0;
		Valid = false;
	}
};

// MW和GF组合观测值数据的定义(用于粗差探测)
struct MWGF_t
{
	short   Prn;
	sys_t   Sys;
	double  MW, GF;
	double  PIF;           //无电离层组合观测值
	int     n;             // 平滑MW组合值的计数器

	MWGF_t()
	{
		Prn = n = 0;
		Sys = sys_t::UNKS;
		MW = GF = PIF = 0.0;
	}
};

// 每颗卫星位置、速度和钟差等的中间计算结果 
struct SatPos_t
{
	double  Pos[3], V[3];        //卫星位置
	double  ClkOft;              //卫星钟差
	double  ClkSft;              //卫星钟速
	double  Tgd1, Tgd2;          //硬件延迟
	double  Elevation, Azimuth; // 高度角，方位角
	double  TropCorr;           // 对流层延迟改正
	bool    Valid;              // false=没有星历或星历过期 true-星历可用
	SatPos_t()
	{
		Pos[0] = Pos[1] = Pos[2] = 0.0;
		V[0] = V[1] = V[2] = 0.0;
		ClkOft = ClkSft = 0.0;
		Elevation = PI / 2.0;
		Azimuth = 0.0;
		Tgd1 = Tgd2 = TropCorr = 0.0;
		Valid = false;
	}
};

//每个历元观测数据
struct epoch_t 
{
	GPSTIME   Time;				   // Current obs time: GPST.
	short     SatNum;			   // 卫星数
	satobs_t  SatObs[MAXCHANNUM];  // 卫星观测数据
	SatPos_t  SatPVT[MAXCHANNUM];  // 卫星位置等计算结果
	MWGF_t    ComObs[MAXCHANNUM];  // 当前历元的组合观测值

	epoch_t()
	{
		SatNum = 0;
	}
};

//每个历元定位结果
struct epochPOS_t
{
	GPSTIME Time;
	BLH     readpos;             //读取的接收机位置bestpos
	double  Pos[3], Vel[3];
	double  PDOP;
	double  SigmaPos, SigmaVel;  //位置，速度的中误差
	int     SatNum;

	epochPOS_t()
	{
		Time.Week = 0;
		Time.SecOfWeek = 0.0;
		readpos = BLH();
		for (int i = 0; i < 3; i++) { Pos[i] = 0.0; Vel[i] = 0.0; }
		PDOP = SigmaPos = SigmaVel = 0.0;
		SatNum = 0;
	}

};

//每颗卫星的单差观测数据定义
struct SDSatOBS_t
{
	short    Prn;
	sys_t    System;
	short    Valid;
	double   dP[2], dL[2];   // 单位为m

	// 站间单差所用基站和流动站的观测数据在其原始观测值数组中的索引
	short    nBas, nRov;

	SDSatOBS_t()
	{
		Prn = nBas = nRov = 0;
		System = sys_t::UNKS;
		dP[0] = dL[0] = dP[1] = dL[1] = 0.0;
		Valid = -1;
	}
};

// 每个历元的单差观测数据定义 
struct SDEpochOBS_t
{
	GPSTIME      Time;
	short        SatNum;
	short        GPSNum;
	short        BDSNum;
	SDSatOBS_t   SdSatObs[MAXCHANNUM];
	MWGF_t       SdCObs[MAXCHANNUM];

	SDEpochOBS_t()
	{
		SatNum = GPSNum = BDSNum = 0;
	}
};

// 双差相关的数据定义 
struct DDOBS_t
{
	int    RefPrn[2], RefPos[2];							      // 参考星卫星号与其在单差数据中的存储位置索引，0=GPS; 1=BDS
	double GPSDDP[2][MAXCHANNUM - 1], BDSDDP[2][MAXCHANNUM - 1];    // GPS和BDS双差伪距观测值数组 (双频)
	double GPSDDL[2][MAXCHANNUM - 1], BDSDDL[2][MAXCHANNUM - 1];    // GPS和BDS双差载波相位观测值数组
	// 另一颗卫星在单差数据中的存储位置索引
	short  otherGPS[MAXCHANNUM - 1], otherBDS[MAXCHANNUM - 1];
	int    SatsNum, DDSatNum[2];							      // 双差观测值数量，0=GPS; 1=BDS
	double dPos[3];											      // 基线向量
	double PDOP;
	bool   RefGPS, RefBDS;									      // 判断GPS和BDS是否找到参考卫星的标签

	DDOBS_t()
	{
		int i;
		PDOP = 0;
		for (i = 0; i < 2; i++)
		{
			DDSatNum[i] = 0;    // 各卫星系统的双差数量
			RefPos[i] = RefPrn[i] = -1;
		}
		SatsNum = 0;              // 双差观测值总数
		dPos[0] = dPos[1] = dPos[2] = 0.0;
		RefGPS = RefBDS = false;
	}
};


// RTK观测数据
struct rtkData_t
{
	epoch_t      basData0;                              //上一时刻基准站数据
	epoch_t      rovData0;                              //上一时刻流动站数据
	epoch_t      basData;                               //本时刻基准站数据
	epoch_t      rovData;                               //本时刻流动站数据
	SDEpochOBS_t SDObs;									//本历元单差观测数据
	DDOBS_t      DDObs;                                 //本历元双差观测数据
	gpseph_t     GPSEph[MAXGPSPRN], BDSEph[MAXBDSPRN];  //GPS和BDS星历
	epochPOS_t   BasPres;                               //基准站参考位置
	epochPOS_t   RovPres;                               //流动站参考位置
	rtkData_t()
	{
		basData0 = epoch_t();
		rovData0 = epoch_t();
		basData = epoch_t();
		rovData = epoch_t();
		SDObs = SDEpochOBS_t();
		DDObs = DDOBS_t();
		GPSEph[MAXGPSPRN] = { gpseph_t() };
		BDSEph[MAXBDSPRN] = { gpseph_t() };
		BasPres = epochPOS_t();
		RovPres = epochPOS_t();
	}
};

// RTK结果数据定义
struct RTKResult_t
{
	/* 浮点解相关 */
	GPSTIME Time;
	int    GPSNum, BDSNum; //卫星数
	double dX_flo[3];      // 基线向量（dx,dy,dz）Rov-Bas
	double dENU_flo[3];    // 基线向量 (dE,dN,dU)
	double RovPos[3];      // 流动站坐标
	double sigma;          // 验后单位权中误差
	double mxyz[3];        // 坐标（基线）分量的中误差
	double mENU[3];        
	double ms;             // 基线长度精度
	double RDOP;
	double RMS;

	/* 模糊度相关 */
	double Ratio;
	int n; // 模糊度维数
	int m; // 模糊度候选解个数，设为2
	double f[2 * MAXCHANNUM];
	double* fa = NULL; // 模糊度浮点解向量 n*1
	double* Qa = NULL; // 模糊度模糊度方差协方差矩阵 n*n
	double* F = NULL; // 模糊度固定解向量，n*m
	double* s = NULL; // 模糊度残差二次型，1*m

	Matrix Qbb = Matrix(3, 3, 0.0);
	Matrix Qaa = Matrix(n, n, 0.0);
	Matrix Qba = Matrix(3, n, 0.0);
	Matrix Qab = Matrix(n, 3, 0.0);
	Matrix a = Matrix(n, 1, 0.0);
	Matrix b = Matrix(3, 1, 0.0);

	Matrix a_fix = Matrix(n, 1, 0.0);
	Matrix b_fix = Matrix(3, 1, 0.0);
	Matrix Qbb_fix = Matrix(3, 3, 0.0);

	double RovPos_fix[3];
	double dX_fix[3];
	double dENU_fix[3];

	int fixed; //是否固定(1---固定成功，0---固定失败)

	/* 精度评定 */
	double dx_FixStd[3];
	double dENU_FixStd[3];

	/* 误差 */
	double delta_x_flo[3];
	double delta_e_flo[3];
	double delta_x_fix[3];
	double delta_e_fix[3];

	RTKResult_t()
	{
		GPSNum = BDSNum = 0;
		dX_flo[0] = dX_flo[1] = dX_flo[2] = 0.0;
		dX_fix[0] = dX_fix[1] = dX_fix[2] = 0.0;
		RovPos_fix[0] = RovPos_fix[1] = RovPos_fix[2] = 0.0;
		mxyz[0] = mxyz[1] = mxyz[2] = 0.0;
		ms = 0.0;
		sigma = 0.0;
		n = 0;
		m = 2;
		for (int i = 0; i < 2 * MAXCHANNUM; i++) f[i] = 0.0;
		Ratio = 0.0;
		fixed = 0;
	}

};

// 配置文件
struct CFGINFO_t
{
	short           datatype;                                      //数据类型（0----文件，1---网口）
	short           systype;									   //系统（0---GPS，1---BDS，2---GPS+BDS）
	double          BasPos_x, BasPos_y, BasPos_z;                  //参考坐标
	char            BasIP[20], RovIP[20];                          //网口IP
	unsigned short  BasPort,RovPort;                               //网口端口
	char            BasInfo[256], RovInfo[256];                    //基站、流动站文件名
	double          dx, dy, dz;                                    //基线向量（xyz）
	double          dE, dN, dU;                                    //基线向量（ENU）

	/*（00-- - GPS单频，01-- - GPS双频，10-- - BDS单频，11-- - BDS双频，20-- - 双系统单频，21-- - 双系统双频）*/
	short           RTKType;                                       //RTK解算类型

	CFGINFO_t()
	{
		BasPos_x = BasPos_y = BasPos_z = 0.0;
		dx = dy = dz = 0.0;
		dE = dN = dU = 0.0;
	}
};

/*
************************************************************
函数名：数据类型转化函数
参数：p 指向所需数据第一个字节的位置
返回值：转化后对应类型的数据
函数功能：将指针p所指向的内存内容解释为对应数据类型的数据
************************************************************
*/
#define U1(p) (*((uint8_t *)(p)))
#define I1(p) (*((int8_t  *)(p)))
static uint16_t U2(uint8_t* p) { uint16_t u; memcpy(&u, p, 2); return u; }
static uint32_t U4(uint8_t* p) { uint32_t u; memcpy(&u, p, 4); return u; }
static int32_t  I4(uint8_t* p) { int32_t  i; memcpy(&i, p, 4); return i; }
static float    R4(uint8_t* p) { float    r; memcpy(&r, p, 4); return r; }
static double   R8(uint8_t* p) { double   r; memcpy(&r, p, 8); return r; }

//解码
//CRC校验函数
uint32_t rtkCRC32(const uint8_t* buff, const int len);
//解码总函数
int DecodeNovOem7Dat(unsigned char Buff[], int& Len, epoch_t* obs, gpseph_t geph[], gpseph_t beph[], epochPOS_t* epos);
//解码Range（观测值）
int decodeRangeb(unsigned char* buff, epoch_t* obs);
//解码GPS星历
int decodeGpsEphem(unsigned char* buff, gpseph_t* geph);
//解码BDS星历
int decodeBdsEphem(unsigned char* buff, gpseph_t* beph);
//解码BestPos
int decodeBestpos(uint8_t* buff, BLH& readpos);

//计算卫星位置
//计算GPS卫星位置、速度、钟差、钟速
bool CalGPSPos(const int Prn, const GPSTIME* t, const gpseph_t* eph, SatPos_t* pos);
//计算BDS卫星位置、速度、钟差、钟速
bool CalBDSPos(const int Prn, const GPSTIME* t, const gpseph_t* eph, SatPos_t* pos);
//计算两种卫星卫星发射时刻的位置、速度、钟差、钟速，并改正位置和速度，计算高度角、方位角、对流层延迟
void CalSatPos(epoch_t* obs, gpseph_t* GPSEph, gpseph_t* BDSEph, XYZ xyz);

//粗差探测
//进行粗差探测并计算IF组合观测值
void DetectError(epoch_t* obs);
//对流层延迟改正（hopfield模型）
double hopfield(double hgt, double elev);

//单点定位
//单点定位SPP
bool SPP(epoch_t* obs, gpseph_t* GPSEph, gpseph_t* BDSEph, epochPOS_t* pos,CFGINFO_t cfg);
//单点测速SPV
void SPV(epoch_t* obs, epochPOS_t* pos);
//计算误差，打印输出结果
void OutputResult(const epochPOS_t* pos, XYZ refpos, Matrix& dENU);

//读取配置文件，获取部分参数设置
bool ReadCFG(const char cfgname[],CFGINFO_t* set);

//rtk
//时间同步（事后）
int TimeSyncPP(ifstream& FBas, ifstream& FRov, rtkData_t* rtkdata);
//时间同步（实时）
int TimeSyncRT(SOCKET& BasSock, SOCKET& RovSock, rtkData_t* rtkdata, vector<epoch_t>& rovOBS, vector<epoch_t>& basOBS, ofstream& BasData, ofstream& RovData);
//存储5个历元时间同步
void timeSyn(vector<epoch_t>rOBS, vector<epoch_t>bOBS, int& irove, int& ibase);
//站间单差函数
void SDEpochObs(epoch_t* RovData, epoch_t* BasData, SDEpochOBS_t* SDObs);
//参考卫星选取函数
void SelectRefSat(const epoch_t* RovData, SDEpochOBS_t* SDObs, DDOBS_t* DDObs);
 //求取双差观测值函数
 void DDEpochObs(SDEpochOBS_t* SDObs, DDOBS_t* DDObs);
 //最小二乘求浮点解函数
 bool LSRTKFloat(rtkData_t* rtkdata, CFGINFO_t cfg, RTKResult_t* result);
 //序贯求固定解函数
 void LSRTKFix(rtkData_t* rtkdata, RTKResult_t* result);
 //输出固定解函数
 void FixOutPut(ofstream& outfile, RTKResult_t* result, CFGINFO_t cfg);