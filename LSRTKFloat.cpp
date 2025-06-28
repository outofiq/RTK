#include "header.h"

/*
****************************************************************************
函数名：最小二乘求浮点解
参数：rtkdata    指向观测数据（双差观测值等）的指针
	  cfg        参数配置信息
	  result     存储解算结果的结构体
返回值：数据同步结果的标志  true-浮点解解算成功  false-浮点解解算失败
函数功能：构建观测方程和随机模型进行最小二乘计算浮点解并进行精度评定；准备模糊度固定需要传入的参数
****************************************************************************
*/
bool LSRTKFloat(rtkData_t* rtkdata, CFGINFO_t cfg, RTKResult_t* result)
{
	result->Time = rtkdata->rovData.Time;
	/* 随机模型常数 */
	double sigma0_2 = 0.003 * 0.003;
	double alpha = 1.5;
	double BasPos[3] = { cfg.BasPos_x,cfg.BasPos_y,cfg.BasPos_z };
	rtkdata->BasPres.Pos[0] = BasPos[0];
	rtkdata->BasPres.Pos[1] = BasPos[1];
	rtkdata->BasPres.Pos[2] = BasPos[2];
	double RovPos[3] = { rtkdata->RovPres.Pos[0],rtkdata->RovPres.Pos[1] ,rtkdata->RovPres.Pos[2] }; //流动站单点定位结果

	//GPS参考卫星位置 (0---基准站，1---流动站）
	double GPSSatPos0[3] = { rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[0]].nBas].Pos[0],
							 rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[0]].nBas].Pos[1],
							 rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[0]].nBas].Pos[2] };
	double GPSSatPos1[3] = { rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[0]].nRov].Pos[0],
							 rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[0]].nRov].Pos[1],
							 rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[0]].nRov].Pos[2] };
	//BDS参考卫星位置 (0---基准站，1---流动站)
	double BDSSatPos0[3] = { rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[1]].nBas].Pos[0],
							 rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[1]].nBas].Pos[1],
							 rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[1]].nBas].Pos[2] };
	double BDSSatPos1[3] = { rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[1]].nRov].Pos[0],
							 rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[1]].nRov].Pos[1],
							 rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[1]].nRov].Pos[2] };

	//基准站到GPS参考卫星的几何距离
	double BasToRefGPS = sqrt((BasPos[0] - GPSSatPos0[0]) * (BasPos[0] - GPSSatPos0[0]) + (BasPos[1] - GPSSatPos0[1]) * (BasPos[1] - GPSSatPos0[1]) + (BasPos[2] - GPSSatPos0[2]) * (BasPos[2] - GPSSatPos0[2]));
	//基准站到BDS参考卫星的几何距离
	double BasToRefBDS = sqrt((BasPos[0] - BDSSatPos0[0]) * (BasPos[0] - BDSSatPos0[0]) + (BasPos[1] - BDSSatPos0[1]) * (BasPos[1] - BDSSatPos0[1]) + (BasPos[2] - BDSSatPos0[2]) * (BasPos[2] - BDSSatPos0[2]));
	//流动站到GPS参考卫星的几何距离
	double RovToRefGPS = sqrt((RovPos[0] - GPSSatPos1[0]) * (RovPos[0] - GPSSatPos1[0]) + (RovPos[1] - GPSSatPos1[1]) * (RovPos[1] - GPSSatPos1[1]) + (RovPos[2] - GPSSatPos1[2]) * (RovPos[2] - GPSSatPos1[2]));
	//流动站到BDS参考卫星的几何距离
	double RovToRefBDS = sqrt((RovPos[0] - BDSSatPos1[0]) * (RovPos[0] - BDSSatPos1[0]) + (RovPos[1] - BDSSatPos1[1]) * (RovPos[1] - BDSSatPos1[1]) + (RovPos[2] - BDSSatPos1[2]) * (RovPos[2] - BDSSatPos1[2]));

	//GPS参考卫星相对于基准站的高度角,观测方差
	double refGPSEB = rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[0]].nBas].Elevation;
	double cosRefGPSEB = cos(refGPSEB);
	//double sinRefGPSEB = sin(refGPSEB);
	double sigma_GPSBas = sigma0_2 * (1 + alpha * cosRefGPSEB * cosRefGPSEB);
	//GPS参考卫星相对于流动站的高度角
	double refGPSER = rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[0]].nRov].Elevation;
	double cosRefGPSER = cos(refGPSER);
	//double sinRefGPSER = sin(refGPSER);
	double sigma_GPSRov = sigma0_2 * (1 + alpha * cosRefGPSER * cosRefGPSER);
	//BDS参考卫星相对于基准站的高度角
	double refBDSEB = rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[1]].nBas].Elevation;
	double cosRefBDSEB = cos(refBDSEB);
	//double sinRefBDSEB = sin(refBDSEB);
	double sigma_BDSBas = sigma0_2 * (1 + alpha * cosRefBDSEB * cosRefBDSEB);
	//BDS参考卫星相对于流动站的高度角
	double refBDSER = rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.RefPos[1]].nRov].Elevation;
	double cosRefBDSER = cos(refBDSER);
	//double sinRefBDSER = sin(refBDSER);
	double sigma_BDSRov = sigma0_2 * (1 + alpha * cosRefBDSER * cosRefBDSER);

	//GPS卫星单差方差
	double GPSDRef = sigma_GPSBas + sigma_GPSRov;
	//BDS卫星单差方差
	double BDSDRef = sigma_BDSBas + sigma_BDSRov;

	double sigma_2 = 0.0; // 单位权中误差
	Matrix Qxyz(3, 3, 0.0); // 待估坐标（基线）的协因数阵
	double m[3] = { 0.0,0.0,0.0 }; // 坐标（基线）分量的中误差
	double dX[3] = { 0.0,0.0,0.0 }; // 基线向量
	Matrix f(1, 3, 0.0);
	double ms = 0.0; // 基线长度精度

	/* GPS单频 */
	if (cfg.RTKType == 00)
	{
		int nGPS = rtkdata->DDObs.DDSatNum[0];
		result->GPSNum = nGPS + 1;
		result->BDSNum = 0;
		int r = 2 * nGPS - 3 - nGPS;
		if (nGPS < 3)
		{
			cerr << "GPS卫星数不足以求解" << endl;
			return false;
		}
		result->n = nGPS;

		Matrix L(2 * nGPS, 1, 0.0);
		Matrix A(2 * nGPS, 3 + nGPS, 0.0);
		Matrix X(3 + nGPS, 1, 0.0);
		Matrix P(2 * nGPS, 2 * nGPS, 0.0);
		Matrix V(2 * nGPS, 1, 0.0);
		Matrix N(3 + nGPS, 3 + nGPS, 0.0);
		Matrix W(3 + nGPS, 1, 0.0);

		double DRef = sigma_GPSBas + sigma_GPSRov; //GPS参考卫星单差方差
		int count = 0;
		bool flag = true; //迭代控制器

		/* 权阵P定义 */
		for (int i = 0; i < nGPS; i++)
		{
			double ER = rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Elevation;
			double EB = rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Elevation;
			double cosER = cos(ER);
			double cosEB = cos(EB);
			double sigmaR2 = sigma0_2 * (1 + alpha * cosER * cosER);
			double sigmaB2 = sigma0_2 * (1 + alpha * cosEB * cosEB);
			double D = sigmaB2 + sigmaR2 + DRef;
			P(i, i) = D;
			P(nGPS + i, nGPS + i) = D * 1e4;
			for (int j = 0; j < nGPS; j++)
			{
				if (i != j) { P(i, j) = DRef; P(nGPS + i, nGPS + j) = DRef * 1e4; }
			}
		}
		P = P.inverse();

		do
		{
			//流动站到GPS参考卫星的几何距离
			double RovToRefGPS = sqrt((RovPos[0] - GPSSatPos1[0]) * (RovPos[0] - GPSSatPos1[0]) + (RovPos[1] - GPSSatPos1[1]) * (RovPos[1] - GPSSatPos1[1]) + (RovPos[2] - GPSSatPos1[2]) * (RovPos[2] - GPSSatPos1[2]));
			
			for (int i = 0; i < nGPS; i++)
			{
				//另一颗GPS卫星在基准站下的坐标
				double BasGPSPos[3] = { rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Pos[0],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Pos[1],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Pos[2] };
				//另一颗GPS卫星在流动站下的坐标
				double RovGPSPos[3] = { rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Pos[0],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Pos[1],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Pos[2] };
				//基准站到另一颗GPS卫星的几何距离
				double BasToGPS = sqrt((BasPos[0] - BasGPSPos[0]) * (BasPos[0] - BasGPSPos[0]) + (BasPos[1] - BasGPSPos[1]) * (BasPos[1] - BasGPSPos[1]) + (BasPos[2] - BasGPSPos[2]) * (BasPos[2] - BasGPSPos[2]));
				//流动站到另一颗GPS卫星的几何距离
				double RovToGPS = sqrt((RovPos[0] - RovGPSPos[0]) * (RovPos[0] - RovGPSPos[0]) + (RovPos[1] - RovGPSPos[1]) * (RovPos[1] - RovGPSPos[1]) + (RovPos[2] - RovGPSPos[2]) * (RovPos[2] - RovGPSPos[2]));
				L(i, 0) = rtkdata->DDObs.GPSDDL[0][i] + RovToRefGPS - BasToRefGPS - RovToGPS + BasToGPS;
				L(nGPS + i, 0) = rtkdata->DDObs.GPSDDP[0][i] + RovToRefGPS - BasToRefGPS - RovToGPS + BasToGPS;

				/* 系数矩阵（A）构建 */
				double aXB = (RovPos[0] - RovGPSPos[0]) / RovToGPS - (RovPos[0] - GPSSatPos1[0]) / RovToRefGPS;
				double aYB = (RovPos[1] - RovGPSPos[1]) / RovToGPS - (RovPos[1] - GPSSatPos1[1]) / RovToRefGPS;
				double aZB = (RovPos[2] - RovGPSPos[2]) / RovToGPS - (RovPos[2] - GPSSatPos1[2]) / RovToRefGPS;
				A(i, 0) = aXB;
				A(i, 1) = aYB;
				A(i, 2) = aZB;
				A(i, i + 3) = WL1_GPS;
				A(nGPS + i, 0) = aXB;
				A(nGPS + i, 1) = aYB;
				A(nGPS + i, 2) = aZB;
			}
			N = A.transpose() * P * A;
			W = A.transpose() * P * L;
			X = N.inverse() * W;
			RovPos[0] = RovPos[0] + X(0, 0);
			RovPos[1] = RovPos[1] + X(1, 0);
			RovPos[2] = RovPos[2] + X(2, 0);

			count++;
			if (fabs(X(0, 0)) < 1e-3 && fabs(X(1, 0)) < 1e-3 && fabs(X(2, 0)) < 1e-3) flag = false;
		} while (flag);
		
		/* 精度评定 */
		V = A * X - L;
		double sigma_2 = (V.transpose() * P * V)(0, 0) / r; // 单位权中误差
		Matrix Qxx = N.inverse();
		Matrix Dxx(3 + result->n, 3 + result->n, 0.0);
		for (int i = 0; i < 3 + result->n; i++)
		{
			for (int j = 0; j < 3 + result->n; j++)
			{
				Dxx(i, j) = sigma_2 * Qxx(i, j);
			}
		}

		result->Qbb = Matrix(3, 3, 0.0);
		result->Qaa = Matrix(result->n, result->n, 0.0);
		result->Qba = Matrix(3, result->n, 0.0);
		result->Qab = Matrix(result->n, 3, 0.0);
		result->a = Matrix(result->n, 1, 0.0);
		result->b = Matrix(3, 1, 0.0);

		result->a_fix = Matrix(result->n, 1, 0.0);
		result->b_fix = Matrix(3, 1, 0.0);
		result->Qbb_fix = Matrix(3, 3, 0.0);
		for (int i = 0; i < 3; i++)
		{
			result->b(i, 0) = X(i, 0);
			for (int j = 0; j < 3; j++)
			{
				Qxyz(i, j) = Qxx(i, j);
				result->Qbb(i, j) = Dxx(i, j);
			}
			for (int k = 0; k < result->n; k++)
			{
				result->Qba(i, k) = Dxx(i, k + 3);
			}
		}

		for (int i = 0; i < result->n; i++)
		{
			result->a(i, 0) = X(i + 3, 0);
			for (int j = 0; j < result->n; j++)
			{
				result->Qaa(i, j) = Dxx(i + 3, j + 3);
			}
			for (int k = 0; k < 3; k++)
			{
				result->Qab(i, k) = Dxx(i + 3, k);
			}
		}
		m[0] = sqrt(sigma_2 * Qxyz(0, 0));
		m[1] = sqrt(sigma_2 * Qxyz(1, 1));
		m[2] = sqrt(sigma_2 * Qxyz(2, 2));
		dX[0] = RovPos[0] - BasPos[0];
		dX[1] = RovPos[1] - BasPos[1];
		dX[2] = RovPos[2] - BasPos[2];
		double S = sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
		f(0, 0) = dX[0] / S;
		f(0, 1) = dX[1] / S;
		f(0, 2) = dX[2] / S;
		double Qs = (f * Qxyz * f.transpose())(0, 0);
		ms = sqrt(sigma_2 * Qs);

		double trQ = 0.0;
		for (int i = 0; i < result->n + 3; i++)
		{
			trQ += Qxx(i, i);
		}
		result->RDOP = sqrt(trQ);

		result->RMS = sqrt((V.transpose() * V)(0, 0) / (2 * nGPS));

		result->RovPos[0] = RovPos[0];
		result->RovPos[1] = RovPos[1];
		result->RovPos[2] = RovPos[2];

		result->dX_flo[0] = dX[0];
		result->dX_flo[1] = dX[1];
		result->dX_flo[2] = dX[2];

		result->sigma = sqrt(sigma_2);

		result->mxyz[0] = m[0];
		result->mxyz[1] = m[1];
		result->mxyz[2] = m[2];

		result->ms = ms;

		/* 模糊度固定所需的参数输入准备 */
		result->fa = mat(result->n, 1);
		result->Qa = mat(result->n, result->n);
		result->F = mat(result->n, result->m);
		result->s = mat(1, result->m);
		for (int i = 0; i < result->n; i++)
		{
			result->fa[i] = X(3 + i, 0);
			for (int j = 0; j < result->n; j++)
			{
				result->Qa[i * result->n + j] = Qxx(3 + i, 3 + j);
			}
		}

		//cout << fixed << setprecision(4) << X(0, 0) << " " << X(1, 0) << " " << X(2, 0) << endl;
		return true;
	}

	/* GPS双频 */
	if (cfg.RTKType == 01)
	{
		int nGPS = rtkdata->DDObs.DDSatNum[0];
		result->GPSNum = nGPS + 1;
		result->BDSNum = 0;
		int r = 4 * nGPS - 3 - 2 * nGPS;
		if (nGPS < 2)
		{
			cerr << "GPS卫星数不足以求解" << endl;
			return false;
		}
		result->n = 2 * nGPS;

		Matrix L(4 * nGPS, 1, 0.0);
		Matrix A(4 * nGPS, 3 + 2 * nGPS, 0.0);
		Matrix X(3 + 2 * nGPS, 1, 0.0);
		Matrix P(4 * nGPS, 4 * nGPS, 0.0);
		Matrix V(4 * nGPS, 1, 0.0);
		Matrix N(3 + 2 * nGPS, 3 + 2 * nGPS, 0.0);
		Matrix W(3 + 2 * nGPS, 1, 0.0);

		double DRef = sigma_GPSBas + sigma_GPSRov; //GPS参考卫星单差方差
		int count = 0;
		bool flag = true; //迭代控制器

		/* 权阵P定义 */
		for (int i = 0; i < nGPS; i++)
		{
			double ER = rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Elevation;
			double EB = rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Elevation;
			double cosER = cos(ER);
			double cosEB = cos(EB);
			double sigmaR2 = sigma0_2 * (1 + alpha * cosER * cosER);
			double sigmaB2 = sigma0_2 * (1 + alpha * cosEB * cosEB);
			double D = sigmaB2 + sigmaR2;
			P(i, i) = D + DRef;
			P(nGPS + i, nGPS + i) = D + DRef;
			P(2 * nGPS + i, 2 * nGPS + i) = (D + DRef) * 1e4;
			P(3 * nGPS + i, 3 * nGPS + i) = (D + DRef) * 1e4;
			for (int j = 0; j < nGPS; j++)
			{
				if (i != j)
				{
					P(i, j) = DRef;
					P(nGPS + i, nGPS + j) = DRef;
					P(2 * nGPS + i, 2 * nGPS + j) = DRef * 1e4;
					P(3 * nGPS + i, 3 * nGPS + j) = DRef * 1e4;
				}
			}
		}
		P = P.inverse();

		do
		{
			//流动站到GPS参考卫星的几何距离
			double RovToRefGPS = sqrt((RovPos[0] - GPSSatPos1[0]) * (RovPos[0] - GPSSatPos1[0]) + (RovPos[1] - GPSSatPos1[1]) * (RovPos[1] - GPSSatPos1[1]) + (RovPos[2] - GPSSatPos1[2]) * (RovPos[2] - GPSSatPos1[2]));
			
			for (int i = 0; i < nGPS; i++)
			{
				//另一颗GPS卫星在基准站下的坐标
				double BasGPSPos[3] = { rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Pos[0],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Pos[1],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Pos[2] };
				//另一颗GPS卫星在流动站下的坐标
				double RovGPSPos[3] = { rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Pos[0],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Pos[1],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Pos[2] };
				//基准站到另一颗GPS卫星的几何距离
				double BasToGPS = sqrt((BasPos[0] - BasGPSPos[0]) * (BasPos[0] - BasGPSPos[0]) + (BasPos[1] - BasGPSPos[1]) * (BasPos[1] - BasGPSPos[1]) + (BasPos[2] - BasGPSPos[2]) * (BasPos[2] - BasGPSPos[2]));
				//流动站到另一颗GPS卫星的几何距离
				double RovToGPS = sqrt((RovPos[0] - RovGPSPos[0]) * (RovPos[0] - RovGPSPos[0]) + (RovPos[1] - RovGPSPos[1]) * (RovPos[1] - RovGPSPos[1]) + (RovPos[2] - RovGPSPos[2]) * (RovPos[2] - RovGPSPos[2]));
				L(i, 0) = rtkdata->DDObs.GPSDDL[0][i] + RovToRefGPS - BasToRefGPS - RovToGPS + BasToGPS;
				L(nGPS + i, 0) = rtkdata->DDObs.GPSDDL[1][i] + RovToRefGPS - BasToRefGPS - RovToGPS + BasToGPS;
				L(2 * nGPS + i, 0) = rtkdata->DDObs.GPSDDP[0][i] + RovToRefGPS - BasToRefGPS - RovToGPS + BasToGPS;
				L(3 * nGPS + i, 0) = rtkdata->DDObs.GPSDDP[1][i] + RovToRefGPS - BasToRefGPS - RovToGPS + BasToGPS;

				/* 系数矩阵（A）构建 */
				double aXB = (RovPos[0] - RovGPSPos[0]) / RovToGPS - (RovPos[0] - GPSSatPos1[0]) / RovToRefGPS;
				double aYB = (RovPos[1] - RovGPSPos[1]) / RovToGPS - (RovPos[1] - GPSSatPos1[1]) / RovToRefGPS;
				double aZB = (RovPos[2] - RovGPSPos[2]) / RovToGPS - (RovPos[2] - GPSSatPos1[2]) / RovToRefGPS;
				A(i, 0) = aXB;
				A(i, 1) = aYB;
				A(i, 2) = aZB;
				A(i, i + 3) = WL1_GPS;
				A(nGPS + i, 0) = aXB;
				A(nGPS + i, 1) = aYB;
				A(nGPS + i, 2) = aZB;
				A(nGPS + i, nGPS + i + 3) = WL2_GPS;
				A(2 * nGPS + i, 0) = aXB;
				A(2 * nGPS + i, 1) = aYB;
				A(2 * nGPS + i, 2) = aZB;
				A(3 * nGPS + i, 0) = aXB;
				A(3 * nGPS + i, 1) = aYB;
				A(3 * nGPS + i, 2) = aZB;
			}
			N = A.transpose() * P * A;
			W = A.transpose() * P * L;
			X = N.inverse() * W;
			RovPos[0] = RovPos[0] + X(0, 0);
			RovPos[1] = RovPos[1] + X(1, 0);
			RovPos[2] = RovPos[2] + X(2, 0);

			count++;
			if (fabs(X(0, 0)) < 1e-3 && fabs(X(1, 0)) < 1e-3 && fabs(X(2, 0)) < 1e-3) flag = false;
		} while (flag);
		
		/* 精度评定 */
		V = A * X - L;
		double sigma_2 = (V.transpose() * P * V)(0, 0) / r; // 单位权中误差
		Matrix Qxx = N.inverse();
		Matrix Dxx(3 + result->n, 3 + result->n, 0.0);
		for (int i = 0; i < 3 + result->n; i++)
		{
			for (int j = 0; j < 3 + result->n; j++)
			{
				Dxx(i, j) = sigma_2 * Qxx(i, j);
			}
		}

		result->Qbb = Matrix(3, 3, 0.0);
		result->Qaa = Matrix(result->n, result->n, 0.0);
		result->Qba = Matrix(3, result->n, 0.0);
		result->Qab = Matrix(result->n, 3, 0.0);
		result->a = Matrix(result->n, 1, 0.0);
		result->b = Matrix(3, 1, 0.0);

		result->a_fix = Matrix(result->n, 1, 0.0);
		result->b_fix = Matrix(3, 1, 0.0);
		result->Qbb_fix = Matrix(3, 3, 0.0);
		for (int i = 0; i < 3; i++)
		{
			result->b(i, 0) = X(i, 0);
			for (int j = 0; j < 3; j++)
			{
				Qxyz(i, j) = Qxx(i, j);
				result->Qbb(i, j) = Dxx(i, j);
			}
			for (int k = 0; k < result->n; k++)
			{
				result->Qba(i, k) = Dxx(i, k + 3);
			}
		}

		for (int i = 0; i < result->n; i++)
		{
			result->a(i, 0) = X(i + 3, 0);
			for (int j = 0; j < result->n; j++)
			{
				result->Qaa(i, j) = Dxx(i + 3, j + 3);
			}
			for (int k = 0; k < 3; k++)
			{
				result->Qab(i, k) = Dxx(i + 3, k);
			}
		}

		m[0] = sqrt(sigma_2 * Qxyz(0, 0));
		m[1] = sqrt(sigma_2 * Qxyz(1, 1));
		m[2] = sqrt(sigma_2 * Qxyz(2, 2));
		dX[0] = RovPos[0] - BasPos[0];
		dX[1] = RovPos[1] - BasPos[1];
		dX[2] = RovPos[2] - BasPos[2];
		double S = sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
		f(0, 0) = dX[0] / S;
		f(0, 1) = dX[1] / S;
		f(0, 2) = dX[2] / S;
		double Qs = (f * Qxyz * f.transpose())(0, 0);
		ms = sqrt(sigma_2 * Qs);

		double trQ = 0.0;
		for (int i = 0; i < result->n + 3; i++)
		{
			trQ += Qxx(i, i);
		}
		result->RDOP = sqrt(trQ);

		result->RMS = sqrt((V.transpose() * V)(0, 0) / (4 * nGPS));

		result->RovPos[0] = RovPos[0];
		result->RovPos[1] = RovPos[1];
		result->RovPos[2] = RovPos[2];

		result->dX_flo[0] = dX[0];
		result->dX_flo[1] = dX[1];
		result->dX_flo[2] = dX[2];

		result->sigma = sqrt(sigma_2);

		result->mxyz[0] = m[0];
		result->mxyz[1] = m[1];
		result->mxyz[2] = m[2];

		result->ms = ms;

		/* 模糊度固定所需的参数输入准备 */
		result->fa = mat(result->n, 1);
		result->Qa = mat(result->n, result->n);
		result->F = mat(result->n, result->m);
		result->s = mat(1, result->m);
		for (int i = 0; i < result->n; i++)
		{
			result->fa[i] = X(3 + i, 0);
			for (int j = 0; j < result->n; j++)
			{
				result->Qa[i * result->n + j] = Qxx(3 + i, 3 + j);
			}
		}

		//cout << fixed << setprecision(4) << X(0, 0) << " " << X(1, 0) << " " << X(2, 0) << endl;
		return true;
	}

	/* BDS单频 */
	if (cfg.RTKType == 10)
	{
		int nBDS = rtkdata->DDObs.DDSatNum[1];
		result->GPSNum = 0;
		result->BDSNum = nBDS + 1;
		int r = 2 * nBDS - 3 - nBDS;
		if (nBDS < 3)
		{
			cerr << "BDS卫星数不足以求解" << endl;
			return false;
		}
		result->n = nBDS;

		Matrix L(2 * nBDS, 1, 0.0);
		Matrix A(2 * nBDS, 3 + nBDS, 0.0);
		Matrix X(3 + nBDS, 1, 0.0);
		Matrix P(2 * nBDS, 2 * nBDS, 0.0);
		Matrix V(2 * nBDS, 1, 0.0);
		Matrix N(3 + nBDS, 3 + nBDS, 0.0);
		Matrix W(3 + nBDS, 1, 0.0);

		double DRef = sigma_BDSBas + sigma_BDSRov; //BDS参考卫星单差方差
		int count = 0;
		bool flag = true; //迭代控制器

		/* 权阵P定义 */
		for (int i = 0; i < nBDS; i++)
		{
			double ER = rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Elevation;
			double EB = rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Elevation;
			double cosER = cos(ER);
			double cosEB = cos(EB);
			double sigmaR2 = sigma0_2 * (1 + alpha * cosER * cosER);
			double sigmaB2 = sigma0_2 * (1 + alpha * cosEB * cosEB);
			double D = sigmaB2 + sigmaR2;
			P(i, i) = D + DRef;
			P(nBDS + i, nBDS + i) = (D + DRef) * 1e4;
			for (int j = 0; j < nBDS; j++)
			{
				if (i != j) { P(i, j) = DRef; P(nBDS + i, nBDS + j) = DRef * 1e4; }
			}
		}
		P = P.inverse();

		do
		{
			//流动站到BDS参考卫星的几何距离
			double RovToRefBDS = sqrt((RovPos[0] - BDSSatPos1[0]) * (RovPos[0] - BDSSatPos1[0]) + (RovPos[1] - BDSSatPos1[1]) * (RovPos[1] - BDSSatPos1[1]) + (RovPos[2] - BDSSatPos1[2]) * (RovPos[2] - BDSSatPos1[2]));
			for (int i = 0; i < nBDS; i++)
			{
				//另一颗BDS卫星在基准站下的坐标
				double BasBDSPos[3] = { rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Pos[0],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Pos[1],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Pos[2] };
				//另一颗BDS卫星在流动站下的坐标
				double RovBDSPos[3] = { rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Pos[0],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Pos[1],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Pos[2] };
				//基准站到另一颗BDS卫星的几何距离
				double BasToBDS = sqrt((BasPos[0] - BasBDSPos[0]) * (BasPos[0] - BasBDSPos[0]) + (BasPos[1] - BasBDSPos[1]) * (BasPos[1] - BasBDSPos[1]) + (BasPos[2] - BasBDSPos[2]) * (BasPos[2] - BasBDSPos[2]));
				//流动站到另一颗BDS卫星的几何距离
				double RovToBDS = sqrt((RovPos[0] - RovBDSPos[0]) * (RovPos[0] - RovBDSPos[0]) + (RovPos[1] - RovBDSPos[1]) * (RovPos[1] - RovBDSPos[1]) + (RovPos[2] - RovBDSPos[2]) * (RovPos[2] - RovBDSPos[2]));
				L(i, 0) = rtkdata->DDObs.BDSDDL[0][i] + RovToRefBDS - BasToRefBDS - RovToBDS + BasToBDS;
				L(nBDS + i, 0) = rtkdata->DDObs.BDSDDP[0][i] + RovToRefBDS - BasToRefBDS - RovToBDS + BasToBDS;

				/* 系数矩阵（A）构建 */
				double aXB = (RovPos[0] - RovBDSPos[0]) / RovToBDS - (RovPos[0] - BDSSatPos1[0]) / RovToRefBDS;
				double aYB = (RovPos[1] - RovBDSPos[1]) / RovToBDS - (RovPos[1] - BDSSatPos1[1]) / RovToRefBDS;
				double aZB = (RovPos[2] - RovBDSPos[2]) / RovToBDS - (RovPos[2] - BDSSatPos1[2]) / RovToRefBDS;
				A(i, 0) = aXB;
				A(i, 1) = aYB;
				A(i, 2) = aZB;
				A(i, i + 3) = WL1_BDS;
				A(nBDS + i, 0) = aXB;
				A(nBDS + i, 1) = aYB;
				A(nBDS + i, 2) = aZB;
			}
			N = A.transpose() * P * A;
			W = A.transpose() * P * L;
			X = N.inverse() * W;
			RovPos[0] = RovPos[0] + X(0, 0);
			RovPos[1] = RovPos[1] + X(1, 0);
			RovPos[2] = RovPos[2] + X(2, 0);

			count++;
			if (fabs(X(0, 0)) < 1e-3 && fabs(X(1, 0)) < 1e-3 && fabs(X(2, 0)) < 1e-3) flag = false;
		} while (flag);
		
		/* 精度评定 */
		V = A * X - L;
		double sigma_2 = (V.transpose() * P * V)(0, 0) / r; // 单位权中误差
		Matrix Qxx = N.inverse();
		Matrix Dxx(3 + result->n, 3 + result->n, 0.0);
		for (int i = 0; i < 3 + result->n; i++)
		{
			for (int j = 0; j < 3 + result->n; j++)
			{
				Dxx(i, j) = sigma_2 * Qxx(i, j);
			}
		}

		result->Qbb = Matrix(3, 3, 0.0);
		result->Qaa = Matrix(result->n, result->n, 0.0);
		result->Qba = Matrix(3, result->n, 0.0);
		result->Qab = Matrix(result->n, 3, 0.0);
		result->a = Matrix(result->n, 1, 0.0);
		result->b = Matrix(3, 1, 0.0);

		result->a_fix = Matrix(result->n, 1, 0.0);
		result->b_fix = Matrix(3, 1, 0.0);
		result->Qbb_fix = Matrix(3, 3, 0.0);
		for (int i = 0; i < 3; i++)
		{
			result->b(i, 0) = X(i, 0);
			for (int j = 0; j < 3; j++)
			{
				Qxyz(i, j) = Qxx(i, j);
				result->Qbb(i, j) = Dxx(i, j);
			}
			for (int k = 0; k < result->n; k++)
			{
				result->Qba(i, k) = Dxx(i, k + 3);
			}
		}

		for (int i = 0; i < result->n; i++)
		{
			result->a(i, 0) = X(i + 3, 0);
			for (int j = 0; j < result->n; j++)
			{
				result->Qaa(i, j) = Dxx(i + 3, j + 3);
			}
			for (int k = 0; k < 3; k++)
			{
				result->Qab(i, k) = Dxx(i + 3, k);
			}
		}

		m[0] = sqrt(sigma_2 * Qxyz(0, 0));
		m[1] = sqrt(sigma_2 * Qxyz(1, 1));
		m[2] = sqrt(sigma_2 * Qxyz(2, 2));
		dX[0] = RovPos[0] - BasPos[0];
		dX[1] = RovPos[1] - BasPos[1];
		dX[2] = RovPos[2] - BasPos[2];
		double S = sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
		f(0, 0) = dX[0] / S;
		f(0, 1) = dX[1] / S;
		f(0, 2) = dX[2] / S;
		double Qs = (f * Qxyz * f.transpose())(0, 0);
		ms = sqrt(sigma_2 * Qs);

		double trQ = 0.0;
		for (int i = 0; i < result->n + 3; i++)
		{
			trQ += Qxx(i, i);
		}
		result->RDOP = sqrt(trQ);

		result->RMS = sqrt((V.transpose() * V)(0, 0) / (2 * nBDS));

		result->RovPos[0] = RovPos[0];
		result->RovPos[1] = RovPos[1];
		result->RovPos[2] = RovPos[2];

		result->dX_flo[0] = dX[0];
		result->dX_flo[1] = dX[1];
		result->dX_flo[2] = dX[2];

		result->sigma = sqrt(sigma_2);

		result->mxyz[0] = m[0];
		result->mxyz[1] = m[1];
		result->mxyz[2] = m[2];

		result->ms = ms;

		/* 模糊度固定所需的参数输入准备 */
		result->fa = mat(result->n, 1);
		result->Qa = mat(result->n, result->n);
		result->F = mat(result->n, result->m);
		result->s = mat(1, result->m);
		for (int i = 0; i < result->n; i++)
		{
			result->fa[i] = X(3 + i, 0);
			for (int j = 0; j < result->n; j++)
			{
				result->Qa[i * result->n + j] = Qxx(3 + i, 3 + j);
			}
		}
		//cout << fixed << setprecision(4) << X(0, 0) << " " << X(1, 0) << " " << X(2, 0) << endl;
		return true;
	}

	/* BDS双频 */
	if (cfg.RTKType == 11)
	{
		int nBDS = rtkdata->DDObs.DDSatNum[1];
		result->GPSNum = 0;
		result->BDSNum = nBDS + 1;
		int r = 4 * nBDS - 3 - 2 * nBDS;
		if (nBDS < 2)
		{
			cerr << "BDS卫星数不足以求解" << endl;
			return false;
		}
		result->n = 2 * nBDS;

		Matrix L(4 * nBDS, 1, 0.0);
		Matrix A(4 * nBDS, 3 + 2 * nBDS, 0.0);
		Matrix X(3 + 2 * nBDS, 1, 0.0);
		Matrix P(4 * nBDS, 4 * nBDS, 0.0);
		Matrix V(4 * nBDS, 1, 0.0);
		Matrix N(3 + 2 * nBDS, 3 + 2 * nBDS, 0.0);
		Matrix W(3 + 2 * nBDS, 1, 0.0);

		double DRef = sigma_BDSBas + sigma_BDSRov; //BDS参考卫星单差方差
		int count = 0;
		bool flag = true; //迭代控制器

		/* 权阵P定义 */
		for (int i = 0; i < nBDS; i++)
		{
			double ER = rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Elevation;
			double EB = rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Elevation;
			double cosER = cos(ER);
			double cosEB = cos(EB);
			double sigmaR2 = sigma0_2 * (1 + alpha * cosER * cosER);
			double sigmaB2 = sigma0_2 * (1 + alpha * cosEB * cosEB);
			double D = sigmaB2 + sigmaR2;
			P(i, i) = D + DRef;
			P(nBDS + i, nBDS + i) = D + DRef;
			P(2 * nBDS + i, 2 * nBDS + i) = (D + DRef) * 1e4;
			P(3 * nBDS + i, 3 * nBDS + i) = (D + DRef) * 1e4;
			for (int j = 0; j < nBDS; j++)
			{
				if (i != j)
				{
					P(i, j) = DRef;
					P(nBDS + i, nBDS + j) = DRef;
					P(2 * nBDS + i, 2 * nBDS + j) = DRef * 1e4;
					P(3 * nBDS + i, 3 * nBDS + j) = DRef * 1e4;
				}
			}
		}
		P = P.inverse();

		do
		{
			//流动站到BDS参考卫星的几何距离
			double RovToRefBDS = sqrt((RovPos[0] - BDSSatPos1[0]) * (RovPos[0] - BDSSatPos1[0]) + (RovPos[1] - BDSSatPos1[1]) * (RovPos[1] - BDSSatPos1[1]) + (RovPos[2] - BDSSatPos1[2]) * (RovPos[2] - BDSSatPos1[2]));
			for (int i = 0; i < nBDS; i++)
			{
				//另一颗BDS卫星在基准站下的坐标
				double BasBDSPos[3] = { rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Pos[0],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Pos[1],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Pos[2] };
				//另一颗BDS卫星在流动站下的坐标
				double RovBDSPos[3] = { rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Pos[0],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Pos[1],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Pos[2] };
				//基准站到另一颗BDS卫星的几何距离
				double BasToBDS = sqrt((BasPos[0] - BasBDSPos[0]) * (BasPos[0] - BasBDSPos[0]) + (BasPos[1] - BasBDSPos[1]) * (BasPos[1] - BasBDSPos[1]) + (BasPos[2] - BasBDSPos[2]) * (BasPos[2] - BasBDSPos[2]));
				//流动站到另一颗BDS卫星的几何距离
				double RovToBDS = sqrt((RovPos[0] - RovBDSPos[0]) * (RovPos[0] - RovBDSPos[0]) + (RovPos[1] - RovBDSPos[1]) * (RovPos[1] - RovBDSPos[1]) + (RovPos[2] - RovBDSPos[2]) * (RovPos[2] - RovBDSPos[2]));
				L(i, 0) = rtkdata->DDObs.BDSDDL[0][i] + RovToRefBDS - BasToRefBDS - RovToBDS + BasToBDS;
				L(nBDS + i, 0) = rtkdata->DDObs.BDSDDL[1][i] + RovToRefBDS - BasToRefBDS - RovToBDS + BasToBDS;
				L(2 * nBDS + i, 0) = rtkdata->DDObs.BDSDDP[0][i] + RovToRefBDS - BasToRefBDS - RovToBDS + BasToBDS;
				L(3 * nBDS + i, 0) = rtkdata->DDObs.BDSDDP[1][i] + RovToRefBDS - BasToRefBDS - RovToBDS + BasToBDS;

				/* 系数矩阵（A）构建 */
				double aXB = (RovPos[0] - RovBDSPos[0]) / RovToBDS - (RovPos[0] - BDSSatPos1[0]) / RovToRefBDS;
				double aYB = (RovPos[1] - RovBDSPos[1]) / RovToBDS - (RovPos[1] - BDSSatPos1[1]) / RovToRefBDS;
				double aZB = (RovPos[2] - RovBDSPos[2]) / RovToBDS - (RovPos[2] - BDSSatPos1[2]) / RovToRefBDS;
				A(i, 0) = aXB;
				A(i, 1) = aYB;
				A(i, 2) = aZB;
				A(i, i + 3) = WL1_BDS;
				A(nBDS + i, 0) = aXB;
				A(nBDS + i, 1) = aYB;
				A(nBDS + i, 2) = aZB;
				A(nBDS + i, nBDS + i + 3) = WL3_BDS;
				A(2 * nBDS + i, 0) = aXB;
				A(2 * nBDS + i, 1) = aYB;
				A(2 * nBDS + i, 2) = aZB;
				A(3 * nBDS + i, 0) = aXB;
				A(3 * nBDS + i, 1) = aYB;
				A(3 * nBDS + i, 2) = aZB;
			}

			N = A.transpose() * P * A;
			W = A.transpose() * P * L;
			X = N.inverse() * W;
			RovPos[0] = RovPos[0] + X(0, 0);
			RovPos[1] = RovPos[1] + X(1, 0);
			RovPos[2] = RovPos[2] + X(2, 0);

			count++;
			if (fabs(X(0, 0)) < 1e-3 && fabs(X(1, 0)) < 1e-3 && fabs(X(2, 0)) < 1e-3) flag = false;
		} while (flag);
		
		/* 精度评定 */
		V = A * X - L;
		double sigma_2 = (V.transpose() * P * V)(0, 0) / r; // 单位权中误差
		Matrix Qxx = N.inverse();
		Matrix Dxx(3 + result->n, 3 + result->n, 0.0);
		for (int i = 0; i < 3 + result->n; i++)
		{
			for (int j = 0; j < 3 + result->n; j++)
			{
				Dxx(i, j) = sigma_2 * Qxx(i, j);
			}
		}

		result->Qbb = Matrix(3, 3, 0.0);
		result->Qaa = Matrix(result->n, result->n, 0.0);
		result->Qba = Matrix(3, result->n, 0.0);
		result->Qab = Matrix(result->n, 3, 0.0);
		result->a = Matrix(result->n, 1, 0.0);
		result->b = Matrix(3, 1, 0.0);

		result->a_fix = Matrix(result->n, 1, 0.0);
		result->b_fix = Matrix(3, 1, 0.0);
		result->Qbb_fix = Matrix(3, 3, 0.0);
		for (int i = 0; i < 3; i++)
		{
			result->b(i, 0) = X(i, 0);
			for (int j = 0; j < 3; j++)
			{
				Qxyz(i, j) = Qxx(i, j);
				result->Qbb(i, j) = Dxx(i, j);
			}
			for (int k = 0; k < result->n; k++)
			{
				result->Qba(i, k) = Dxx(i, k + 3);
			}
		}

		for (int i = 0; i < result->n; i++)
		{
			result->a(i, 0) = X(i + 3, 0);
			for (int j = 0; j < result->n; j++)
			{
				result->Qaa(i, j) = Dxx(i + 3, j + 3);
			}
			for (int k = 0; k < 3; k++)
			{
				result->Qab(i, k) = Dxx(i + 3, k);
			}
		}

		m[0] = sqrt(sigma_2 * Qxyz(0, 0));
		m[1] = sqrt(sigma_2 * Qxyz(1, 1));
		m[2] = sqrt(sigma_2 * Qxyz(2, 2));
		dX[0] = RovPos[0] - BasPos[0];
		dX[1] = RovPos[1] - BasPos[1];
		dX[2] = RovPos[2] - BasPos[2];
		double S = sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
		f(0, 0) = dX[0] / S;
		f(0, 1) = dX[1] / S;
		f(0, 2) = dX[2] / S;
		double Qs = (f * Qxyz * f.transpose())(0, 0);
		ms = sqrt(sigma_2 * Qs);

		double trQ = 0.0;
		for (int i = 0; i < result->n + 3; i++)
		{
			trQ += Qxx(i, i);
		}
		result->RDOP = sqrt(trQ);

		result->RMS = sqrt((V.transpose() * V)(0, 0) / (4 * nBDS));

		result->RovPos[0] = RovPos[0];
		result->RovPos[1] = RovPos[1];
		result->RovPos[2] = RovPos[2];

		result->dX_flo[0] = dX[0];
		result->dX_flo[1] = dX[1];
		result->dX_flo[2] = dX[2];

		result->sigma = sqrt(sigma_2);

		result->mxyz[0] = m[0];
		result->mxyz[1] = m[1];
		result->mxyz[2] = m[2];

		result->ms = ms;

		/* 模糊度固定所需的参数输入准备 */
		result->fa = mat(result->n, 1);
		result->Qa = mat(result->n, result->n);
		result->F = mat(result->n, result->m);
		result->s = mat(1, result->m);
		for (int i = 0; i < result->n; i++)
		{
			result->fa[i] = X(3 + i, 0);
			for (int j = 0; j < result->n; j++)
			{
				result->Qa[i * result->n + j] = Qxx(3 + i, 3 + j);
			}
		}
		//cout << fixed << setprecision(4) << X(0, 0) << " " << X(1, 0) << " " << X(2, 0) << endl;
		return true;
	}

	/* 双系统单频 */
	if (cfg.RTKType == 20)
	{
		int nGPS = rtkdata->DDObs.DDSatNum[0];
		int nBDS = rtkdata->DDObs.DDSatNum[1];
		result->GPSNum = nGPS + 1;
		result->BDSNum = nBDS + 1;
		int n = rtkdata->DDObs.SatsNum;
		int r = 2 * n - 3 - n;
		if (n < 3)
		{
			cerr << "卫星数不足以求解" << endl;
			return false;
		}
		result->n = n;

		Matrix L(2 * n, 1, 0.0);
		Matrix A(2 * n, 3 + n, 0.0);
		Matrix X(3 + n, 1, 0.0);
		Matrix P(2 * n, 2 * n, 0.0);
		Matrix V(2 * n, 1, 0.0);
		Matrix N(3 + n, 3 + n, 0.0);
		Matrix W(3 + n, 1, 0.0);

		double GPSDRef = sigma_GPSBas + sigma_GPSRov; //GPS参考卫星单差方差
		double BDSDRef = sigma_BDSBas + sigma_BDSRov; //BDS参考卫星单差方差
		int count = 0;
		bool flag = true; //迭代控制器

		/* 权阵P定义 */
		for (int i = 0; i < nGPS; i++)
		{
			double ER = rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Elevation;
			double EB = rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Elevation;
			double cosER = cos(ER);
			double cosEB = cos(EB);
			double sigmaR2 = sigma0_2 * (1 + alpha * cosER * cosER);
			double sigmaB2 = sigma0_2 * (1 + alpha * cosEB * cosEB);
			double D = sigmaB2 + sigmaR2;
			P(i, i) = D + GPSDRef;
			P(n + i, n + i) = (D + GPSDRef) * 1e4;
			for (int j = 0; j < nGPS; j++)
			{
				if (i != j)
				{
					P(i, j) = GPSDRef;
					P(n + i, n + j) = GPSDRef * 1e4;
				}
			}
		}
		for (int i = 0; i < nBDS; i++)
		{
			double ER = rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Elevation;
			double EB = rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Elevation;
			double cosER = cos(ER);
			double cosEB = cos(EB);
			double sigmaR2 = sigma0_2 * (1 + alpha * cosER * cosER);
			double sigmaB2 = sigma0_2 * (1 + alpha * cosEB * cosEB);
			double D = sigmaB2 + sigmaR2;
			P(nGPS + i, nGPS + i) = D + BDSDRef;
			P(n + nGPS + i, n + nGPS + i) = (D + BDSDRef) * 1e4;
			for (int j = 0; j < nBDS; j++)
			{
				if (i != j)
				{
					P(nGPS + i, nGPS + j) = BDSDRef;
					P(n + nGPS + i, n + nGPS + j) = BDSDRef * 1e4;
				}
			}
		}
		P = P.inverse();

		do
		{
			//流动站到GPS参考卫星的几何距离
			double RovToRefGPS = sqrt((RovPos[0] - GPSSatPos1[0]) * (RovPos[0] - GPSSatPos1[0]) + (RovPos[1] - GPSSatPos1[1]) * (RovPos[1] - GPSSatPos1[1]) + (RovPos[2] - GPSSatPos1[2]) * (RovPos[2] - GPSSatPos1[2]));
			//流动站到BDS参考卫星的几何距离
			double RovToRefBDS = sqrt((RovPos[0] - BDSSatPos1[0]) * (RovPos[0] - BDSSatPos1[0]) + (RovPos[1] - BDSSatPos1[1]) * (RovPos[1] - BDSSatPos1[1]) + (RovPos[2] - BDSSatPos1[2]) * (RovPos[2] - BDSSatPos1[2]));
			for (int i = 0; i < nGPS; i++)
			{
				//另一颗GPS卫星在基准站下的坐标
				double BasGPSPos[3] = { rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Pos[0],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Pos[1],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Pos[2] };
				//另一颗GPS卫星在流动站下的坐标
				double RovGPSPos[3] = { rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Pos[0],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Pos[1],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Pos[2] };
				//基准站到另一颗GPS卫星的几何距离
				double BasToGPS = sqrt((BasPos[0] - BasGPSPos[0]) * (BasPos[0] - BasGPSPos[0]) + (BasPos[1] - BasGPSPos[1]) * (BasPos[1] - BasGPSPos[1]) + (BasPos[2] - BasGPSPos[2]) * (BasPos[2] - BasGPSPos[2]));
				//流动站到另一颗GPS卫星的几何距离
				double RovToGPS = sqrt((RovPos[0] - RovGPSPos[0]) * (RovPos[0] - RovGPSPos[0]) + (RovPos[1] - RovGPSPos[1]) * (RovPos[1] - RovGPSPos[1]) + (RovPos[2] - RovGPSPos[2]) * (RovPos[2] - RovGPSPos[2]));
				L(i, 0) = rtkdata->DDObs.GPSDDL[0][i] + RovToRefGPS - BasToRefGPS - RovToGPS + BasToGPS;
				L(n + i, 0) = rtkdata->DDObs.GPSDDP[0][i] + RovToRefGPS - BasToRefGPS - RovToGPS + BasToGPS;

				/* 系数矩阵（A）构建 */
				double aXB = (RovPos[0] - RovGPSPos[0]) / RovToGPS - (RovPos[0] - GPSSatPos1[0]) / RovToRefGPS;
				double aYB = (RovPos[1] - RovGPSPos[1]) / RovToGPS - (RovPos[1] - GPSSatPos1[1]) / RovToRefGPS;
				double aZB = (RovPos[2] - RovGPSPos[2]) / RovToGPS - (RovPos[2] - GPSSatPos1[2]) / RovToRefGPS;
				A(i, 0) = aXB;
				A(i, 1) = aYB;
				A(i, 2) = aZB;
				A(i, i + 3) = WL1_GPS;
				A(n + i, 0) = aXB;
				A(n + i, 1) = aYB;
				A(n + i, 2) = aZB;
			}

			for (int i = 0; i < nBDS; i++)
			{
				//另一颗BDS卫星在基准站下的坐标
				double BasBDSPos[3] = { rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Pos[0],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Pos[1],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Pos[2] };
				//另一颗BDS卫星在流动站下的坐标
				double RovBDSPos[3] = { rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Pos[0],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Pos[1],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Pos[2] };
				//基准站到另一颗BDS卫星的几何距离
				double BasToBDS = sqrt((BasPos[0] - BasBDSPos[0]) * (BasPos[0] - BasBDSPos[0]) + (BasPos[1] - BasBDSPos[1]) * (BasPos[1] - BasBDSPos[1]) + (BasPos[2] - BasBDSPos[2]) * (BasPos[2] - BasBDSPos[2]));
				//流动站到另一颗BDS卫星的几何距离
				double RovToBDS = sqrt((RovPos[0] - RovBDSPos[0]) * (RovPos[0] - RovBDSPos[0]) + (RovPos[1] - RovBDSPos[1]) * (RovPos[1] - RovBDSPos[1]) + (RovPos[2] - RovBDSPos[2]) * (RovPos[2] - RovBDSPos[2]));
				L(nGPS + i, 0) = rtkdata->DDObs.BDSDDL[0][i] + RovToRefBDS - BasToRefBDS - RovToBDS + BasToBDS;
				L(n + nGPS + i, 0) = rtkdata->DDObs.BDSDDP[0][i] + RovToRefBDS - BasToRefBDS - RovToBDS + BasToBDS;

				/* 系数矩阵（A）构建 */
				double aXB = (RovPos[0] - RovBDSPos[0]) / RovToBDS - (RovPos[0] - BDSSatPos1[0]) / RovToRefBDS;
				double aYB = (RovPos[1] - RovBDSPos[1]) / RovToBDS - (RovPos[1] - BDSSatPos1[1]) / RovToRefBDS;
				double aZB = (RovPos[2] - RovBDSPos[2]) / RovToBDS - (RovPos[2] - BDSSatPos1[2]) / RovToRefBDS;
				A(nGPS + i, 0) = aXB;
				A(nGPS + i, 1) = aYB;
				A(nGPS + i, 2) = aZB;
				A(nGPS + i, nGPS + i + 3) = WL1_BDS;
				A(n + nGPS + i, 0) = aXB;
				A(n + nGPS + i, 1) = aYB;
				A(n + nGPS + i, 2) = aZB;
			}
			N = A.transpose() * P * A;
			W = A.transpose() * P * L;
			X = N.inverse() * W;
			RovPos[0] = RovPos[0] + X(0, 0);
			RovPos[1] = RovPos[1] + X(1, 0);
			RovPos[2] = RovPos[2] + X(2, 0);

			count++;
			if (fabs(X(0, 0)) < 1e-3 && fabs(X(1, 0)) < 1e-3 && fabs(X(2, 0)) < 1e-3) flag = false;
		} while (flag);

		/* 精度评定 */
		V = A * X - L;
		double sigma_2 = (V.transpose() * P * V)(0, 0) / r; // 单位权中误差
		Matrix Qxx = N.inverse();
		Matrix Dxx(3 + result->n, 3 + result->n, 0.0);
		for (int i = 0; i < 3 + result->n; i++)
		{
			for (int j = 0; j < 3 + result->n; j++)
			{
				Dxx(i, j) = sigma_2 * Qxx(i, j);
			}
		}

		result->Qbb = Matrix(3, 3, 0.0);
		result->Qaa = Matrix(result->n, result->n, 0.0);
		result->Qba = Matrix(3, result->n, 0.0);
		result->Qab = Matrix(result->n, 3, 0.0);
		result->a = Matrix(result->n, 1, 0.0);
		result->b = Matrix(3, 1, 0.0);

		result->a_fix = Matrix(result->n, 1, 0.0);
		result->b_fix = Matrix(3, 1, 0.0);
		result->Qbb_fix = Matrix(3, 3, 0.0);
		for (int i = 0; i < 3; i++)
		{
			result->b(i, 0) = X(i, 0);
			for (int j = 0; j < 3; j++)
			{
				Qxyz(i, j) = Qxx(i, j);
				result->Qbb(i, j) = Dxx(i, j);
			}
			for (int k = 0; k < result->n; k++)
			{
				result->Qba(i, k) = Dxx(i, k + 3);
			}
		}

		for (int i = 0; i < result->n; i++)
		{
			result->a(i, 0) = X(i + 3, 0);
			for (int j = 0; j < result->n; j++)
			{
				result->Qaa(i, j) = Dxx(i + 3, j + 3);
			}
			for (int k = 0; k < 3; k++)
			{
				result->Qab(i, k) = Dxx(i + 3, k);
			}
		}

		m[0] = sqrt(sigma_2 * Qxyz(0, 0));
		m[1] = sqrt(sigma_2 * Qxyz(1, 1));
		m[2] = sqrt(sigma_2 * Qxyz(2, 2));
		dX[0] = RovPos[0] - BasPos[0];
		dX[1] = RovPos[1] - BasPos[1];
		dX[2] = RovPos[2] - BasPos[2];
		double S = sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
		f(0, 0) = dX[0] / S;
		f(0, 1) = dX[1] / S;
		f(0, 2) = dX[2] / S;
		double Qs = (f * Qxyz * f.transpose())(0, 0);
		ms = sqrt(sigma_2 * Qs);

		double trQ = 0.0;
		for (int i = 0; i < result->n + 3; i++)
		{
			trQ += Qxx(i, i);
		}
		result->RDOP = sqrt(trQ);

		result->RMS = sqrt((V.transpose() * V)(0, 0) / (2 * n));

		result->RovPos[0] = RovPos[0];
		result->RovPos[1] = RovPos[1];
		result->RovPos[2] = RovPos[2];

		result->dX_flo[0] = dX[0];
		result->dX_flo[1] = dX[1];
		result->dX_flo[2] = dX[2];

		result->sigma = sqrt(sigma_2);

		result->mxyz[0] = m[0];
		result->mxyz[1] = m[1];
		result->mxyz[2] = m[2];

		result->ms = ms;
		
		/* 模糊度固定所需的参数输入准备 */
		result->fa = mat(result->n, 1);
		result->Qa = mat(result->n, result->n);
		result->F = mat(result->n, result->m);
		result->s = mat(1, result->m);
		for (int i = 0; i < result->n; i++)
		{
			result->fa[i] = X(3 + i, 0);
			for (int j = 0; j < result->n; j++)
			{
				result->Qa[i * result->n + j] = Qxx(3 + i, 3 + j);
			}
		}
		//cout << fixed << setprecision(4) << X(0, 0) << " " << X(1, 0) << " " << X(2, 0) << endl;
		return true;
	}

	/* 双系统双频 */
	if (cfg.RTKType == 21)
	{
		int nGPS = rtkdata->DDObs.DDSatNum[0];
		int nBDS = rtkdata->DDObs.DDSatNum[1];
		result->GPSNum = nGPS + 1;
		result->BDSNum = nBDS + 1;
		int n = rtkdata->DDObs.SatsNum;
		int r = 4 * n - 3 - 2 * n;
		if (n < 2)
		{
			cerr << "卫星数不足以求解" << endl;
			return false;
		}
		result->n = 2 * n;

		Matrix L(4 * n, 1, 0.0);
		Matrix V(4 * n, 1, 0.0);
		Matrix A(4 * n, 3 + 2 * n, 0.0);
		Matrix X(3 + 2 * n, 1, 0.0);
		Matrix P(4 * n, 4 * n, 0.0);
		Matrix N(3 + 2 * n, 3 + 2 * n, 0.0);
		Matrix W(3 + 2 * n, 1, 0.0);

		int count = 0; //迭代计数器
		bool flag = true; //迭代控制器


		/* 权阵P的定义（由于每次迭代不改变，故放在循环外面） */
		for (int i = 0; i < nGPS; i++)
		{
			double GPSER = rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Elevation;
			double GPSEB = rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Elevation;
			double cosGPSER = cos(GPSER);
			double cosGPSEB = cos(GPSEB);
			double sigmaR2 = sigma0_2 * (1 + alpha * cosGPSER * cosGPSER);
			double sigmaB2 = sigma0_2 * (1 + alpha * cosGPSEB * cosGPSEB);
			double DGPS = sigmaB2 + sigmaR2 + GPSDRef;
			P(i, i) = DGPS;
			P(nGPS + i, nGPS + i) = DGPS;
			P(2 * n + i, 2 * n + i) = DGPS * 1e4;
			P(2 * n + nGPS + i, 2 * n + nGPS + i) = DGPS * 1e4;
			for (int j = 0; j < nGPS; j++)
			{
				if (i != j)
				{
					P(i, j) = GPSDRef;
					P(nGPS + i, nGPS + j) = GPSDRef;
					P(2 * n + i, 2 * n + j) = GPSDRef * 1e4;
					P(2 * n + nGPS + i, 2 * n + nGPS + j) = GPSDRef * 1e4;
				}
			}
		}
		for (int i = 0; i < nBDS; i++)
		{
			double BDSER = rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Elevation;
			double BDSEB = rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Elevation;
			double cosBDSER = cos(BDSER);
			double cosBDSEB = cos(BDSEB);
			double sigmaR2 = sigma0_2 * (1 + alpha * cosBDSER * cosBDSER);
			double sigmaB2 = sigma0_2 * (1 + alpha * cosBDSEB * cosBDSEB);
			double DBDS = sigmaB2 + sigmaR2 + BDSDRef;
			P(2 * nGPS + i, 2 * nGPS + i) = DBDS;
			P(2 * nGPS + nBDS + i, 2 * nGPS + nBDS + i) = DBDS;
			P(2 * n + 2 * nGPS + i, 2 * n + 2 * nGPS + i) = DBDS * 1e4;
			P(2 * n + 2 * nGPS + nBDS + i, 2 * n + 2 * nGPS + nBDS + i) = DBDS * 1e4;
			for (int j = 0; j < nBDS; j++)
			{
				if (i != j)
				{
					P(2 * nGPS + i, 2 * nGPS + j) = BDSDRef;
					P(2 * nGPS + nBDS + i, 2 * nGPS + nBDS + j) = BDSDRef;
					P(2 * n + 2 * nGPS + i, 2 * n + 2 * nGPS + j) = BDSDRef * 1e4;
					P(2 * n + 2 * nGPS + nBDS + i, 2 * n + 2 * nGPS + nBDS + j) = BDSDRef * 1e4;
				}
			}
		}
		P = P.inverse();

		do
		{
			//流动站到GPS参考卫星的几何距离
			double RovToRefGPS = sqrt((RovPos[0] - GPSSatPos1[0]) * (RovPos[0] - GPSSatPos1[0]) + (RovPos[1] - GPSSatPos1[1]) * (RovPos[1] - GPSSatPos1[1]) + (RovPos[2] - GPSSatPos1[2]) * (RovPos[2] - GPSSatPos1[2]));
			//流动站到BDS参考卫星的几何距离
			double RovToRefBDS = sqrt((RovPos[0] - BDSSatPos1[0]) * (RovPos[0] - BDSSatPos1[0]) + (RovPos[1] - BDSSatPos1[1]) * (RovPos[1] - BDSSatPos1[1]) + (RovPos[2] - BDSSatPos1[2]) * (RovPos[2] - BDSSatPos1[2]));
			for (int i = 0; i < nGPS; i++)
			{

				//另一颗GPS卫星在基准站下的坐标
				double BasGPSPos[3] = { rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Pos[0],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Pos[1],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nBas].Pos[2] };
				//另一颗GPS卫星在流动站下的坐标
				double RovGPSPos[3] = { rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Pos[0],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Pos[1],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherGPS[i]].nRov].Pos[2] };
				//基准站到另一颗GPS卫星的几何距离
				double BasToGPS = sqrt((BasPos[0] - BasGPSPos[0]) * (BasPos[0] - BasGPSPos[0]) + (BasPos[1] - BasGPSPos[1]) * (BasPos[1] - BasGPSPos[1]) + (BasPos[2] - BasGPSPos[2]) * (BasPos[2] - BasGPSPos[2]));
				//流动站到另一颗GPS卫星的几何距离
				double RovToGPS = sqrt((RovPos[0] - RovGPSPos[0]) * (RovPos[0] - RovGPSPos[0]) + (RovPos[1] - RovGPSPos[1]) * (RovPos[1] - RovGPSPos[1]) + (RovPos[2] - RovGPSPos[2]) * (RovPos[2] - RovGPSPos[2]));
				L(i, 0) = rtkdata->DDObs.GPSDDL[0][i] + RovToRefGPS - BasToRefGPS - RovToGPS + BasToGPS;
				L(nGPS + i, 0) = rtkdata->DDObs.GPSDDL[1][i] + RovToRefGPS - BasToRefGPS - RovToGPS + BasToGPS;
				L(2 * n + i, 0) = rtkdata->DDObs.GPSDDP[0][i] + RovToRefGPS - BasToRefGPS - RovToGPS + BasToGPS;
				L(2 * n + nGPS + i, 0) = rtkdata->DDObs.GPSDDP[1][i] + RovToRefGPS - BasToRefGPS - RovToGPS + BasToGPS;

				/* 系数矩阵（A）构建 */
				double aXBGPS = (RovPos[0] - RovGPSPos[0]) / RovToGPS - (RovPos[0] - GPSSatPos1[0]) / RovToRefGPS;
				double aYBGPS = (RovPos[1] - RovGPSPos[1]) / RovToGPS - (RovPos[1] - GPSSatPos1[1]) / RovToRefGPS;
				double aZBGPS = (RovPos[2] - RovGPSPos[2]) / RovToGPS - (RovPos[2] - GPSSatPos1[2]) / RovToRefGPS;
				A(i, 0) = aXBGPS;
				A(i, 1) = aYBGPS;
				A(i, 2) = aZBGPS;
				A(i, i + 3) = WL1_GPS;
				A(nGPS + i, 0) = aXBGPS;
				A(nGPS + i, 1) = aYBGPS;
				A(nGPS + i, 2) = aZBGPS;
				A(nGPS + i, nGPS + i + 3) = WL2_GPS;
				A(2 * n + i, 0) = aXBGPS;
				A(2 * n + i, 1) = aYBGPS;
				A(2 * n + i, 2) = aZBGPS;
				A(2 * n + nGPS + i, 0) = aXBGPS;
				A(2 * n + nGPS + i, 1) = aYBGPS;
				A(2 * n + nGPS + i, 2) = aZBGPS;
			}

			for (int i = 0; i < nBDS; i++)
			{
				//另一颗BDS卫星在基准站下的坐标
				double BasBDSPos[3] = { rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Pos[0],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Pos[1],
										rtkdata->basData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nBas].Pos[2] };
				//另一颗BDS卫星在流动站下的坐标
				double RovBDSPos[3] = { rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Pos[0],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Pos[1],
										rtkdata->rovData.SatPVT[rtkdata->SDObs.SdSatObs[rtkdata->DDObs.otherBDS[i]].nRov].Pos[2] };
				//基准站到另一颗BDS卫星的几何距离
				double BasToBDS = sqrt((BasPos[0] - BasBDSPos[0]) * (BasPos[0] - BasBDSPos[0]) + (BasPos[1] - BasBDSPos[1]) * (BasPos[1] - BasBDSPos[1]) + (BasPos[2] - BasBDSPos[2]) * (BasPos[2] - BasBDSPos[2]));
				//流动站到另一颗BDS卫星的几何距离
				double RovToBDS = sqrt((RovPos[0] - RovBDSPos[0]) * (RovPos[0] - RovBDSPos[0]) + (RovPos[1] - RovBDSPos[1]) * (RovPos[1] - RovBDSPos[1]) + (RovPos[2] - RovBDSPos[2]) * (RovPos[2] - RovBDSPos[2]));
				L(2 * nGPS + i, 0) = rtkdata->DDObs.BDSDDL[0][i] + RovToRefBDS - BasToRefBDS - RovToBDS + BasToBDS;
				L(2 * nGPS + nBDS + i, 0) = rtkdata->DDObs.BDSDDL[1][i] + RovToRefBDS - BasToRefBDS - RovToBDS + BasToBDS;
				L(2 * n + 2 * nGPS + i, 0) = rtkdata->DDObs.BDSDDP[0][i] + RovToRefBDS - BasToRefBDS - RovToBDS + BasToBDS;
				L(2 * n + 2 * nGPS + nBDS + i, 0) = rtkdata->DDObs.BDSDDP[1][i] + RovToRefBDS - BasToRefBDS - RovToBDS + BasToBDS;

				/* 系数矩阵（A）构建 */
				double aXBBDS = (RovPos[0] - RovBDSPos[0]) / RovToBDS - (RovPos[0] - BDSSatPos1[0]) / RovToRefBDS;
				double aYBBDS = (RovPos[1] - RovBDSPos[1]) / RovToBDS - (RovPos[1] - BDSSatPos1[1]) / RovToRefBDS;
				double aZBBDS = (RovPos[2] - RovBDSPos[2]) / RovToBDS - (RovPos[2] - BDSSatPos1[2]) / RovToRefBDS;
				A(2 * nGPS + i, 0) = aXBBDS;
				A(2 * nGPS + i, 1) = aYBBDS;
				A(2 * nGPS + i, 2) = aZBBDS;
				A(2 * nGPS + i, 2 * nGPS + i + 3) = WL1_BDS;
				A(2 * nGPS + nBDS + i, 0) = aXBBDS;
				A(2 * nGPS + nBDS + i, 1) = aYBBDS;
				A(2 * nGPS + nBDS + i, 2) = aZBBDS;
				A(2 * nGPS + nBDS + i, 2 * nGPS + nBDS + i + 3) = WL3_BDS;
				A(2 * n + 2 * nGPS + i, 0) = aXBBDS;
				A(2 * n + 2 * nGPS + i, 1) = aYBBDS;
				A(2 * n + 2 * nGPS + i, 2) = aZBBDS;
				A(2 * n + 2 * nGPS + nBDS + i, 0) = aXBBDS;
				A(2 * n + 2 * nGPS + nBDS + i, 1) = aYBBDS;
				A(2 * n + 2 * nGPS + nBDS + i, 2) = aZBBDS;
			}
			Matrix A_T = A.transpose();
			N = A_T * P * A;
			Matrix Q = N.inverse();
			W = A_T * P * L;
			X = Q * W;
			RovPos[0] = RovPos[0] + X(0, 0);
			RovPos[1] = RovPos[1] + X(1, 0);
			RovPos[2] = RovPos[2] + X(2, 0);

			count++;
			if (fabs(X(0, 0)) < 0.001 && fabs(X(1, 0)) < 0.001 && fabs(X(2, 0)) < 0.001) flag = false;
		} while (flag);

		/* 精度评定 */
		V = A * X - L;
		double sigma_2 = (V.transpose() * P * V)(0, 0) / r; // 单位权中误差
		Matrix Qxx = N.inverse();
		Matrix Dxx(3 + result->n, 3 + result->n, 0.0);
		for (int i = 0; i < 3 + result->n; i++)
		{
			for (int j = 0; j < 3 + result->n; j++)
			{
				Dxx(i, j) = sigma_2 * Qxx(i, j);
			}
		}

		result->Qbb = Matrix(3, 3, 0.0);
		result->Qaa = Matrix(result->n, result->n, 0.0);
		result->Qba = Matrix(3, result->n, 0.0);
		result->Qab = Matrix(result->n, 3, 0.0);
		result->a = Matrix(result->n, 1, 0.0);
		result->b = Matrix(3, 1, 0.0);

		result->a_fix = Matrix(result->n, 1, 0.0);
		result->b_fix = Matrix(3, 1, 0.0);
		result->Qbb_fix = Matrix(3, 3, 0.0);
		for (int i = 0; i < 3; i++)
		{
			result->b(i, 0) = X(i, 0);
			for (int j = 0; j < 3; j++)
			{
				Qxyz(i, j) = Qxx(i, j);
				result->Qbb(i, j) = Dxx(i, j);
			}
			for (int k = 0; k < result->n; k++)
			{
				result->Qba(i, k) = Dxx(i, k + 3);
			}
		}

		for (int i = 0; i < result->n; i++)
		{
			result->a(i, 0) = X(i + 3, 0);
			for (int j = 0; j < result->n; j++)
			{
				result->Qaa(i, j) = Dxx(i + 3, j + 3);
			}
			for (int k = 0; k < 3; k++)
			{
				result->Qab(i, k) = Dxx(i + 3, k);
			}
		}
		m[0] = sqrt(sigma_2 * Qxyz(0, 0));
		m[1] = sqrt(sigma_2 * Qxyz(1, 1));
		m[2] = sqrt(sigma_2 * Qxyz(2, 2));
		dX[0] = RovPos[0] - BasPos[0];
		dX[1] = RovPos[1] - BasPos[1];
		dX[2] = RovPos[2] - BasPos[2];
		double S = sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
		f(0, 0) = dX[0] / S;
		f(0, 1) = dX[1] / S;
		f(0, 2) = dX[2] / S;
		double Qs = (f * Qxyz * f.transpose())(0, 0);
		ms = sqrt(sigma_2 * Qs);

		double trQ = 0.0;
		for (int i = 0; i < result->n + 3; i++)
		{
			trQ += Qxx(i, i);
		}
		result->RDOP = sqrt(trQ);

		result->RMS = sqrt((V.transpose() * V)(0, 0) / (4 * n));

		result->RovPos[0] = RovPos[0];
		result->RovPos[1] = RovPos[1];
		result->RovPos[2] = RovPos[2];

		result->dX_flo[0] = dX[0];
		result->dX_flo[1] = dX[1];
		result->dX_flo[2] = dX[2];

		result->sigma = sqrt(sigma_2);

		result->mxyz[0] = m[0];
		result->mxyz[1] = m[1];
		result->mxyz[2] = m[2];

		result->ms = ms;

		/* 模糊度固定所需的参数输入准备 */
		result->fa = mat(result->n, 1);
		result->Qa = mat(result->n, result->n);
		result->F = mat(result->n, result->m);
		result->s = mat(1, result->m);
		for (int i = 0; i < result->n; i++)
		{
			result->fa[i] = X(3 + i, 0);
			for (int j = 0; j < result->n; j++)
			{
				result->Qa[i * result->n + j] = Qxx(3 + i, 3 + j);
			}
		}
		
		//std::cout << fixed << setprecision(4) << RovPos[0] << " " << RovPos[1] << " " << RovPos[2];
		//std::cout << fixed << setprecision(4) << result->sigma << " " << m[0] << " " << m[1] << " " << m[2] << " " << ms << endl;
		return true;
	}
}