#include "timeTrans.h"

/*
**********************************************************
函数名：取小数部分
参数：Value    一个浮点数
返回值：Value的小数部分
函数功能：获得一个浮点数的小数部分
**********************************************************
*/
double FRAC(const double Value)
{
	return (Value - (int)Value);
}

/*
**********************************************************
函数名：通用时转儒略日
参数：Common    通用时
      MJD       儒略日
返回值：true—转成功
函数功能：将通用时转化为儒略日
**********************************************************
*/
bool CommonToMJD(const COMMONTIME Common, MJDTIME* MJD)
{
	int y, m;
	double JD, UT, TempMJD;
	if (Common.Month <= 2)
	{
		y = Common.Year - 1;
		m = Common.Month + 12;
	}
	if (Common.Month > 2)
	{
		y = Common.Year;
		m = Common.Month;
	}
	UT = Common.Hour + Common.Minute / 60.0 + Common.Second / 3600.0;
	JD = (int)(365.25 * y) + (int)(30.6001 * (m + 1)) + Common.Day + UT / 24.0 + 1720981.5;
	TempMJD = JD - 2400000.5;
	MJD->Days = (int)(TempMJD);
	MJD->FracDay = (TempMJD)-(MJD->Days);
	return true;
}

/*
**********************************************************
函数名：儒略日转化为通用时
参数：MJD     儒略日
      Common  通用时
返回值：true—转化成功
函数功能：将儒略日转化为通用时
**********************************************************
*/
bool MJDToCommon(const MJDTIME MJD, COMMONTIME* Common)
{
	int a, b, c, d, e;
	double UTDay, UT;
	a = (int)(MJD.Days + MJD.FracDay + 2400000.5 + 0.5);
	b = a + 1537;
	c = (int)((b - 122.1) / 365.25);
	d = (int)(365.25 * c);
	e = (int)((b - d) * 1.0 / 30.6001);
	UTDay = b - d - (int)(30.6001 * e) + FRAC(MJD.Days + MJD.FracDay + 2400000.5 + 0.5);
	Common->Day = (int)(UTDay);
	Common->Month = e - 1 - 12 * (int)(e * 1.0 / 14.0);
	Common->Year = c - 4715 - (int)((7 + Common->Month) * 1.0 / 10.0);
	UT = (FRAC(UTDay) * 24.0);
	Common->Hour = int(MJD.FracDay * 24);
	Common->Minute = int((MJD.FracDay * 24 - Common->Hour) * 60);
	Common->Second = MJD.FracDay * 24 * 3600 - Common->Hour * 3600 - Common->Minute * 60;
	return true;
}

/*
**********************************************************
函数名：儒略日转GPS时
参数：MJD    儒略日
      GPS    GPS时
返回值：true—转化成功
函数功能：将儒略日转化为GPS时
**********************************************************
*/
bool MJDToGPS(const MJDTIME MJD, GPSTIME* GPS)
{
	GPS->Week = (int)((MJD.Days + MJD.FracDay - 44244) / 7.0);
	GPS->SecOfWeek = (MJD.Days + MJD.FracDay - 44244 - GPS->Week * 7) * 86400.0;
	return true;
}

/*
**********************************************************
函数名：GPS时转儒略日
参数：GPS    GPS时
      MJD    儒略日
返回值：true—转化成功
函数功能：将GPS时转化为儒略日
**********************************************************
*/
bool GPSToMJD(const GPSTIME GPS, MJDTIME* MJD)
{
	double TempMJD;
	TempMJD = 44244 + GPS.Week * 7 + GPS.SecOfWeek / 86400.0;
	MJD->Days = (int)(TempMJD);
	MJD->FracDay = (TempMJD)-(MJD->Days);
	return true;
}

/*
**********************************************************
函数名：求解两个GPS时的时间差
参数：T1    第一个GPS时
	  T2    第二个GPS时
返回值：时间差（以s为单位）
函数功能：求解两个GPS时的时间差，以秒为单位输出
**********************************************************
*/
double diffTime(const GPSTIME* T1, const GPSTIME* T2)
{
	double diffTime = (T1->Week - T2->Week) * 7 * 24 * 3600 + T1->SecOfWeek - T2->SecOfWeek;
	return diffTime;
}
