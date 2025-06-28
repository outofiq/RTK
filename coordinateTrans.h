#include <cmath>

// �������
#define R_WGS84  6378137.0                // WGS84����ĳ����� m 
#define F_WGS84  1.0/298.257223563        // WGS84����ı���
#define Omega_WGS 7.2921151467e-5         // WGS84�������ת���ٶ� rad/s
#define GM_WGS   398600.5e+9              // WGS84���������������GM m3/s2

#define R_CGS2K  6378137.0                // CGCS2000����ĳ����� m
#define F_CGS2K  1.0/298.257222101        // CGCS2000����ı���
#define Omega_BDS 7.2921150e-5            // CGCS2000�������ת���ٶ� rad/s
#define GM_BDS   398600.4418e+9           // CGCS2000���������������GM m3/s2

	struct XYZ
	{
		double x;
		double y;
		double z;
		XYZ()
		{
			x = 0.0;
			y = 0.0;
			z = 0.0;
		}
	};

	struct BLH
	{
		double B;
		double L;
		double H;
		BLH()
		{
			B = 0.0;
			L = 0.0;
			H = 0.0;
		}
	};

	//XYZת��ΪBLH
	void XYZToBLH(const XYZ Descartes, BLH* Geodetic, const double R, const double F);
	//BLHת��ΪXYZ
	void BLHToXYZ(const BLH Geodetic, XYZ* Descartes, const double R, const double F);