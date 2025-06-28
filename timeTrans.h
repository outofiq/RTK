#include <cmath>

//ͨ��ʱ�䶨��
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

	//�������ն���
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

	//GPSʱ�䶨��
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

	//ͨ��ʱת��Ϊ��������
	bool CommonToMJD(const COMMONTIME Common, MJDTIME* MJD);
	//��������ת��Ϊͨ��ʱ
	bool MJDToCommon(const MJDTIME MJD, COMMONTIME* Common);
	//��������ת��ΪGPSʱ
	bool MJDToGPS(const MJDTIME MJD, GPSTIME* GPS);
	//GPSʱת��Ϊ��������
	bool GPSToMJD(const GPSTIME GPS, MJDTIME* MJD);
	//ȡС�����ֺ���
	double FRAC(const double Value);
	//��GPSʱ��s��
	double diffTime(const GPSTIME* T1, const GPSTIME* T2);