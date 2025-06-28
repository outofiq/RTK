#include "coordinateTrans.h"

/*
**********************************************************
函数名：笛卡尔系转大地坐标系
参数：Descartes    输入的笛卡尔坐标XYZ
	  Geodetic	   转化得到的大地坐标BLH
	  R            椭球长半轴
	  F            椭球扁率
函数功能：将笛卡尔坐标XYZ转化为大地坐标BLH
**********************************************************
*/
void XYZToBLH(const XYZ Descartes, BLH* Geodetic,const double R,const double F)
{
	double x = Descartes.x;
	double y = Descartes.y;
	double z = Descartes.z;

	const double e2 = 2.0 * F - F * F;

	double S = sqrt(x * x + y * y);
	double L = atan2(y, x);
	double B = 0;
	double N = 0;
	double tempB = atan2(z, S);

	int counter = 0;
	while (fabs(B - tempB) > 1e-12 && counter < 25)
	{
		B = tempB;
		N = R / sqrt(1 - e2 * sin(B) * sin(B));
		tempB = atan2(z + N * e2 * sin(B), S);
		counter++;
	}

	Geodetic->B = B;
	Geodetic->L = L;
	Geodetic->H = S / cos(B) - N;
}

/*
**********************************************************
函数名：大地坐标系转笛卡尔坐标系
参数：Geodetic     输入的大地坐标BLH
	  Descartes	   转化得到的笛卡尔坐标XYZ
	  R            椭球长半轴
	  F            椭球扁率
函数功能：将大地坐标坐标BLH转化为笛卡尔坐标XYZ
**********************************************************
*/
void BLHToXYZ(const BLH Geodetic, XYZ* Descartes,const double R,const double F)
{
	double B = Geodetic.B;
	double L = Geodetic.L;
	double H = Geodetic.H;
	double e2 = 2 * F - F * F;
	double N = R / sqrt(1 - e2 * sin(B) * sin(B));
	Descartes->x = (N + H) * cos(B) * cos(L);
	Descartes->y = (N + H) * cos(B) * sin(L);
	Descartes->z = (N * (1 - e2) + H) * sin(B);
}