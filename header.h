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
#define MAXGEOPRN       5            /* ����GEO���Ǻ� */ //�����ƣ�BDS-3��GEO���ų�������
#define MAXBDSPRN       63
#define MAXRAWLEN       40960
#define MAXBUFF	        20480

/* math common uints */
#define PI          3.1415926535898   /* PI */  
#define D2R			(PI/180.0)        /* radians to degree */
#define R2D         (180.0/PI)        /* degree to radians */
#define CLIGHT		299792458.0       /* Speed of light  [m/s]; IAU 1976  */

// GPS�����źŵ�һЩ����
#define  FG1_GPS  1575.42E6           // L1�ź�Ƶ�� 
#define  FG2_GPS  1227.60E6           // L2�ź�Ƶ��
#define  WL1_GPS  (CLIGHT/FG1_GPS)    // L1�źŲ���
#define  WL2_GPS  (CLIGHT/FG2_GPS)    // L2�źŲ���

// BDS�����źŵ�һЩ����
#define  FG1_BDS  1561.098E6          // B1�źŵ�Ƶ�� 
#define  FG3_BDS  1268.520E6          // B3�źŵ�Ƶ�� 
#define  WL1_BDS  (CLIGHT/FG1_BDS)    // B1�źŲ���
#define  WL3_BDS  (CLIGHT/FG3_BDS)    // B3�źŲ���

#define SECPERHOUR 3600.0             /* Seconds per hour */
#define SECPERDAY  86400.0            /* Seconds per day */
#define SECPERWEEK 604800.0           /* Seconds per week */

/* ����ϵͳ���� ---------------------------------------------------------------*/
enum class sys_t { UNKS = 0, GPS, BDS };

//GPS+BDS�㲥����
struct gpseph_t 
{
	unsigned short PRN;
	sys_t          Sys;
	GPSTIME        TOC, TOE; //s
	short	       health; //0-������1-������
	double	       ClkBias, ClkDrift, ClkDriftRate; //�Ӳ� ���� ��Ʈ
	double	       IODE, IODC; //������������ ��������������
	double         TGD1, TGD2; //�㲥����������ʱ��Ⱥ�ӳ�
	double	       SqrtA, ecc, M0, OMEGA0, i0, omega, OMEGADot, iDot, DeltaN;
	double	       Crs, Cuc, Cus, Cic, Cis, Crc; //�㶯��
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

//ÿ�����ǹ۲�����
struct satobs_t 
{
	sys_t    Sys;
	short    Prn;
	//GPS��L1-0��L2-1   BDS��B1I-0��B3I-1
	double   p[2],l[2],d[2]; //α�� ��λ ������
	double   CN0[2];         //����� 
	bool     Valid;          // ���壨��˫Ƶ�����ڣ���Ч��


	satobs_t()
	{
		Prn = 0;
		Sys = sys_t::UNKS;
		p[0] = p[1] = l[0] = l[1] = d[0] = d[1] = 0.0;
		CN0[0] = CN0[1] = 0.0;
		Valid = false;
	}
};

// MW��GF��Ϲ۲�ֵ���ݵĶ���(���ڴֲ�̽��)
struct MWGF_t
{
	short   Prn;
	sys_t   Sys;
	double  MW, GF;
	double  PIF;           //�޵������Ϲ۲�ֵ
	int     n;             // ƽ��MW���ֵ�ļ�����

	MWGF_t()
	{
		Prn = n = 0;
		Sys = sys_t::UNKS;
		MW = GF = PIF = 0.0;
	}
};

// ÿ������λ�á��ٶȺ��Ӳ�ȵ��м������ 
struct SatPos_t
{
	double  Pos[3], V[3];        //����λ��
	double  ClkOft;              //�����Ӳ�
	double  ClkSft;              //��������
	double  Tgd1, Tgd2;          //Ӳ���ӳ�
	double  Elevation, Azimuth; // �߶Ƚǣ���λ��
	double  TropCorr;           // �������ӳٸ���
	bool    Valid;              // false=û���������������� true-��������
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

//ÿ����Ԫ�۲�����
struct epoch_t 
{
	GPSTIME   Time;				   // Current obs time: GPST.
	short     SatNum;			   // ������
	satobs_t  SatObs[MAXCHANNUM];  // ���ǹ۲�����
	SatPos_t  SatPVT[MAXCHANNUM];  // ����λ�õȼ�����
	MWGF_t    ComObs[MAXCHANNUM];  // ��ǰ��Ԫ����Ϲ۲�ֵ

	epoch_t()
	{
		SatNum = 0;
	}
};

//ÿ����Ԫ��λ���
struct epochPOS_t
{
	GPSTIME Time;
	BLH     readpos;             //��ȡ�Ľ��ջ�λ��bestpos
	double  Pos[3], Vel[3];
	double  PDOP;
	double  SigmaPos, SigmaVel;  //λ�ã��ٶȵ������
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

//ÿ�����ǵĵ���۲����ݶ���
struct SDSatOBS_t
{
	short    Prn;
	sys_t    System;
	short    Valid;
	double   dP[2], dL[2];   // ��λΪm

	// վ�䵥�����û�վ������վ�Ĺ۲���������ԭʼ�۲�ֵ�����е�����
	short    nBas, nRov;

	SDSatOBS_t()
	{
		Prn = nBas = nRov = 0;
		System = sys_t::UNKS;
		dP[0] = dL[0] = dP[1] = dL[1] = 0.0;
		Valid = -1;
	}
};

// ÿ����Ԫ�ĵ���۲����ݶ��� 
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

// ˫����ص����ݶ��� 
struct DDOBS_t
{
	int    RefPrn[2], RefPos[2];							      // �ο������Ǻ������ڵ��������еĴ洢λ��������0=GPS; 1=BDS
	double GPSDDP[2][MAXCHANNUM - 1], BDSDDP[2][MAXCHANNUM - 1];    // GPS��BDS˫��α��۲�ֵ���� (˫Ƶ)
	double GPSDDL[2][MAXCHANNUM - 1], BDSDDL[2][MAXCHANNUM - 1];    // GPS��BDS˫���ز���λ�۲�ֵ����
	// ��һ�������ڵ��������еĴ洢λ������
	short  otherGPS[MAXCHANNUM - 1], otherBDS[MAXCHANNUM - 1];
	int    SatsNum, DDSatNum[2];							      // ˫��۲�ֵ������0=GPS; 1=BDS
	double dPos[3];											      // ��������
	double PDOP;
	bool   RefGPS, RefBDS;									      // �ж�GPS��BDS�Ƿ��ҵ��ο����ǵı�ǩ

	DDOBS_t()
	{
		int i;
		PDOP = 0;
		for (i = 0; i < 2; i++)
		{
			DDSatNum[i] = 0;    // ������ϵͳ��˫������
			RefPos[i] = RefPrn[i] = -1;
		}
		SatsNum = 0;              // ˫��۲�ֵ����
		dPos[0] = dPos[1] = dPos[2] = 0.0;
		RefGPS = RefBDS = false;
	}
};


// RTK�۲�����
struct rtkData_t
{
	epoch_t      basData0;                              //��һʱ�̻�׼վ����
	epoch_t      rovData0;                              //��һʱ������վ����
	epoch_t      basData;                               //��ʱ�̻�׼վ����
	epoch_t      rovData;                               //��ʱ������վ����
	SDEpochOBS_t SDObs;									//����Ԫ����۲�����
	DDOBS_t      DDObs;                                 //����Ԫ˫��۲�����
	gpseph_t     GPSEph[MAXGPSPRN], BDSEph[MAXBDSPRN];  //GPS��BDS����
	epochPOS_t   BasPres;                               //��׼վ�ο�λ��
	epochPOS_t   RovPres;                               //����վ�ο�λ��
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

// RTK������ݶ���
struct RTKResult_t
{
	/* �������� */
	GPSTIME Time;
	int    GPSNum, BDSNum; //������
	double dX_flo[3];      // ����������dx,dy,dz��Rov-Bas
	double dENU_flo[3];    // �������� (dE,dN,dU)
	double RovPos[3];      // ����վ����
	double sigma;          // ���λȨ�����
	double mxyz[3];        // ���꣨���ߣ������������
	double mENU[3];        
	double ms;             // ���߳��Ⱦ���
	double RDOP;
	double RMS;

	/* ģ������� */
	double Ratio;
	int n; // ģ����ά��
	int m; // ģ���Ⱥ�ѡ���������Ϊ2
	double f[2 * MAXCHANNUM];
	double* fa = NULL; // ģ���ȸ�������� n*1
	double* Qa = NULL; // ģ����ģ���ȷ���Э������� n*n
	double* F = NULL; // ģ���ȹ̶���������n*m
	double* s = NULL; // ģ���Ȳв�����ͣ�1*m

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

	int fixed; //�Ƿ�̶�(1---�̶��ɹ���0---�̶�ʧ��)

	/* �������� */
	double dx_FixStd[3];
	double dENU_FixStd[3];

	/* ��� */
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

// �����ļ�
struct CFGINFO_t
{
	short           datatype;                                      //�������ͣ�0----�ļ���1---���ڣ�
	short           systype;									   //ϵͳ��0---GPS��1---BDS��2---GPS+BDS��
	double          BasPos_x, BasPos_y, BasPos_z;                  //�ο�����
	char            BasIP[20], RovIP[20];                          //����IP
	unsigned short  BasPort,RovPort;                               //���ڶ˿�
	char            BasInfo[256], RovInfo[256];                    //��վ������վ�ļ���
	double          dx, dy, dz;                                    //����������xyz��
	double          dE, dN, dU;                                    //����������ENU��

	/*��00-- - GPS��Ƶ��01-- - GPS˫Ƶ��10-- - BDS��Ƶ��11-- - BDS˫Ƶ��20-- - ˫ϵͳ��Ƶ��21-- - ˫ϵͳ˫Ƶ��*/
	short           RTKType;                                       //RTK��������

	CFGINFO_t()
	{
		BasPos_x = BasPos_y = BasPos_z = 0.0;
		dx = dy = dz = 0.0;
		dE = dN = dU = 0.0;
	}
};

/*
************************************************************
����������������ת������
������p ָ���������ݵ�һ���ֽڵ�λ��
����ֵ��ת�����Ӧ���͵�����
�������ܣ���ָ��p��ָ����ڴ����ݽ���Ϊ��Ӧ�������͵�����
************************************************************
*/
#define U1(p) (*((uint8_t *)(p)))
#define I1(p) (*((int8_t  *)(p)))
static uint16_t U2(uint8_t* p) { uint16_t u; memcpy(&u, p, 2); return u; }
static uint32_t U4(uint8_t* p) { uint32_t u; memcpy(&u, p, 4); return u; }
static int32_t  I4(uint8_t* p) { int32_t  i; memcpy(&i, p, 4); return i; }
static float    R4(uint8_t* p) { float    r; memcpy(&r, p, 4); return r; }
static double   R8(uint8_t* p) { double   r; memcpy(&r, p, 8); return r; }

//����
//CRCУ�麯��
uint32_t rtkCRC32(const uint8_t* buff, const int len);
//�����ܺ���
int DecodeNovOem7Dat(unsigned char Buff[], int& Len, epoch_t* obs, gpseph_t geph[], gpseph_t beph[], epochPOS_t* epos);
//����Range���۲�ֵ��
int decodeRangeb(unsigned char* buff, epoch_t* obs);
//����GPS����
int decodeGpsEphem(unsigned char* buff, gpseph_t* geph);
//����BDS����
int decodeBdsEphem(unsigned char* buff, gpseph_t* beph);
//����BestPos
int decodeBestpos(uint8_t* buff, BLH& readpos);

//��������λ��
//����GPS����λ�á��ٶȡ��Ӳ����
bool CalGPSPos(const int Prn, const GPSTIME* t, const gpseph_t* eph, SatPos_t* pos);
//����BDS����λ�á��ٶȡ��Ӳ����
bool CalBDSPos(const int Prn, const GPSTIME* t, const gpseph_t* eph, SatPos_t* pos);
//���������������Ƿ���ʱ�̵�λ�á��ٶȡ��Ӳ���٣�������λ�ú��ٶȣ�����߶Ƚǡ���λ�ǡ��������ӳ�
void CalSatPos(epoch_t* obs, gpseph_t* GPSEph, gpseph_t* BDSEph, XYZ xyz);

//�ֲ�̽��
//���дֲ�̽�Ⲣ����IF��Ϲ۲�ֵ
void DetectError(epoch_t* obs);
//�������ӳٸ�����hopfieldģ�ͣ�
double hopfield(double hgt, double elev);

//���㶨λ
//���㶨λSPP
bool SPP(epoch_t* obs, gpseph_t* GPSEph, gpseph_t* BDSEph, epochPOS_t* pos,CFGINFO_t cfg);
//�������SPV
void SPV(epoch_t* obs, epochPOS_t* pos);
//��������ӡ������
void OutputResult(const epochPOS_t* pos, XYZ refpos, Matrix& dENU);

//��ȡ�����ļ�����ȡ���ֲ�������
bool ReadCFG(const char cfgname[],CFGINFO_t* set);

//rtk
//ʱ��ͬ�����º�
int TimeSyncPP(ifstream& FBas, ifstream& FRov, rtkData_t* rtkdata);
//ʱ��ͬ����ʵʱ��
int TimeSyncRT(SOCKET& BasSock, SOCKET& RovSock, rtkData_t* rtkdata, vector<epoch_t>& rovOBS, vector<epoch_t>& basOBS, ofstream& BasData, ofstream& RovData);
//�洢5����Ԫʱ��ͬ��
void timeSyn(vector<epoch_t>rOBS, vector<epoch_t>bOBS, int& irove, int& ibase);
//վ�䵥���
void SDEpochObs(epoch_t* RovData, epoch_t* BasData, SDEpochOBS_t* SDObs);
//�ο�����ѡȡ����
void SelectRefSat(const epoch_t* RovData, SDEpochOBS_t* SDObs, DDOBS_t* DDObs);
 //��ȡ˫��۲�ֵ����
 void DDEpochObs(SDEpochOBS_t* SDObs, DDOBS_t* DDObs);
 //��С�����󸡵�⺯��
 bool LSRTKFloat(rtkData_t* rtkdata, CFGINFO_t cfg, RTKResult_t* result);
 //�����̶��⺯��
 void LSRTKFix(rtkData_t* rtkdata, RTKResult_t* result);
 //����̶��⺯��
 void FixOutPut(ofstream& outfile, RTKResult_t* result, CFGINFO_t cfg);