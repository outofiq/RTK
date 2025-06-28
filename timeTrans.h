#include <cmath>

//通用时间定义
	struct COMMONTIME
	{
		unsigned short Year;
		unsigned short Month;
		unsigned short Day;
		unsigned short Hour;
		unsigned short Minute;
		double Second;

		COMMONTIME()
		{
			Year = 0;
			Month = 0;
			Day = 0;
			Hour = 0;
			Minute = 0;
			Second = 0.0;
		}

	};

	//简化儒略日定义
	struct MJDTIME
	{
		int Days;
		double FracDay;

		MJDTIME()
		{
			Days = 0;
			FracDay = 0.0;
		}
	};

	//GPS时间定义
	struct GPSTIME
	{
		unsigned short Week;
		double SecOfWeek;

		GPSTIME()
		{
			Week = 0;
			SecOfWeek = 0.0;
		}
	};

	//通用时转化为简化儒略日
	bool CommonToMJD(const COMMONTIME Common, MJDTIME* MJD);
	//简化儒略日转换为通用时
	bool MJDToCommon(const MJDTIME MJD, COMMONTIME* Common);
	//简化儒略日转化为GPS时
	bool MJDToGPS(const MJDTIME MJD, GPSTIME* GPS);
	//GPS时转化为简化儒略日
	bool GPSToMJD(const GPSTIME GPS, MJDTIME* MJD);
	//取小数部分函数
	double FRAC(const double Value);
	//求GPS时间差（s）
	double diffTime(const GPSTIME* T1, const GPSTIME* T2);