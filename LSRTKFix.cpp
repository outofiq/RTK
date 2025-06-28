#include "header.h"

/*
****************************************************************************
函数名：RTK固定解更新
参数：rtkdata    指向数据存储处的指针（两个站的观测值、星历、参考定位结果）
	  result     存储RTK解算结果的结构体
函数功能：模糊度固定后更新RTK解算结果，得到的解为固定解
****************************************************************************
*/
void LSRTKFix(rtkData_t* rtkdata,RTKResult_t* result)
{
	Matrix b = result->b;
	Matrix Qbb = result->Qbb;
	Matrix Qaa = result->Qaa;
	Matrix Qba = result->Qba;
	Matrix Qab = result->Qab;
	Matrix a = result->a;
	Matrix a_fix = result->a_fix;
	Matrix Qaa_1 = Qaa.inverse();

	if (result->Ratio>2.0)
	{
		result->b_fix = b - Qba * Qaa_1 * (a - a_fix);
		result->Qbb_fix = Qbb - Qba * Qaa_1 * Qab;
	}
	else if(result->Ratio<=2.0)
	{
		cout << "固定失败" << endl;
		result->b_fix = b;
		result->Qbb_fix = Qbb;
	}
	result->RovPos_fix[0] = result->RovPos[0] + result->b_fix(0, 0);
	result->RovPos_fix[1] = result->RovPos[1] + result->b_fix(1, 0);
	result->RovPos_fix[2] = result->RovPos[2] + result->b_fix(2, 0);

	result->dX_fix[0] = result->RovPos_fix[0] - rtkdata->BasPres.Pos[0];
	result->dX_fix[1] = result->RovPos_fix[1] - rtkdata->BasPres.Pos[1];
	result->dX_fix[2] = result->RovPos_fix[2] - rtkdata->BasPres.Pos[2];

	result->dx_FixStd[0] = sqrt(result->Qbb_fix(0, 0));
	result->dx_FixStd[1] = sqrt(result->Qbb_fix(1, 1));
	result->dx_FixStd[2] = sqrt(result->Qbb_fix(2, 2));
}

void FixOutPut(ofstream& outfile, RTKResult_t* result, CFGINFO_t cfg)
{
	XYZ xyz;
	xyz.x = cfg.BasPos_x;
	xyz.y = cfg.BasPos_y;
	xyz.z = cfg.BasPos_z;
	BLH blh;
	XYZToBLH(xyz, &blh, R_CGS2K, F_CGS2K);
	Matrix Mat(3, 3, 0);
	Mat(0, 0) = -sin(blh.L);
	Mat(0, 1) = cos(blh.L);
	Mat(0, 2) = 0;
	Mat(1, 0) = -sin(blh.B) * cos(blh.L);
	Mat(1, 1) = -sin(blh.B) * sin(blh.L);
	Mat(1, 2) = cos(blh.B);
	Mat(2, 0) = cos(blh.B) * cos(blh.L);
	Mat(2, 1) = cos(blh.B) * sin(blh.L);
	Mat(2, 2) = sin(blh.B);

	Matrix dx_flo(3, 1, 0);
	dx_flo(0, 0) = result->dX_flo[0];
	dx_flo(1, 0) = result->dX_flo[1];
	dx_flo(2, 0) = result->dX_flo[2];
	Matrix dENU_flo = Mat * dx_flo;
	result->dENU_flo[0] = dENU_flo(0, 0);
	result->dENU_flo[1] = dENU_flo(1, 0);
	result->dENU_flo[2] = dENU_flo(2, 0);

	Matrix mxyz(3, 1, 0);
	mxyz(0, 0) = result->mxyz[0];
	mxyz(1, 0) = result->mxyz[1];
	mxyz(2, 0) = result->mxyz[2];
	Matrix mENU = Mat * mxyz;
	result->mENU[0] = mENU(0, 0);
	result->mENU[1] = mENU(1, 0);
	result->mENU[2] = mENU(2, 0);

	Matrix dX_fix(3, 1, 0);
	dX_fix(0, 0) = result->dX_fix[0];
	dX_fix(1, 0) = result->dX_fix[1];
	dX_fix(2, 0) = result->dX_fix[2];
	Matrix dENU_fix = Mat * dX_fix;
	result->dENU_fix[0] = dENU_fix(0, 0);
	result->dENU_fix[1] = dENU_fix(1, 0);
	result->dENU_fix[2] = dENU_fix(2, 0);

	Matrix dX_FixStd(3, 1, 0);
	dX_FixStd(0, 0) = result->dx_FixStd[0];
	dX_FixStd(1, 0) = result->dx_FixStd[1];
	dX_FixStd(2, 0) = result->dx_FixStd[2];
	Matrix dENU_FixStd = Mat * dX_FixStd;
	result->dENU_FixStd[0] = dENU_FixStd(0, 0);
	result->dENU_FixStd[1] = dENU_FixStd(1, 0);
	result->dENU_FixStd[2] = dENU_FixStd(2, 0);

	result->delta_x_flo[0] = result->dX_flo[0] - cfg.dx;
	result->delta_x_flo[1] = result->dX_flo[1] - cfg.dy;
	result->delta_x_flo[2] = result->dX_flo[2] - cfg.dz;

	result->delta_e_flo[0] = result->dENU_flo[0] - cfg.dE;
	result->delta_e_flo[1] = result->dENU_flo[1] - cfg.dN;
	result->delta_e_flo[2] = result->dENU_flo[2] - cfg.dU;

	result->delta_x_fix[0] = result->dX_fix[0] - cfg.dx;
	result->delta_x_fix[1] = result->dX_fix[1] - cfg.dy;
	result->delta_x_fix[2] = result->dX_fix[2] - cfg.dz;

	result->delta_e_fix[0] = result->dENU_fix[0] - cfg.dE;
	result->delta_e_fix[1] = result->dENU_fix[1] - cfg.dN;
	result->delta_e_fix[2] = result->dENU_fix[2] - cfg.dU;

	if (result->Ratio>2.0)
	{
		cout << result->Time.Week << " " << fixed << setprecision(3) << result->Time.SecOfWeek << " Fix " << "GPSNum=" << result->GPSNum << " BDSNum=" << result->BDSNum << endl;
		outfile << result->Time.Week << " " << fixed << setprecision(3) << result->Time.SecOfWeek << " Fix " << " dX(m)= " <<
			fixed << setprecision(4) <<setw(8)<< result->dX_fix[0] << "  dY(m)= " <<setw(8)<< result->dX_fix[1] << "  dZ(m)= "<<setw(8) << result->dX_fix[2] <<
			" delta_X_fix= " << setw(8) <<result->delta_x_fix[0]<< " delta_Y_fix= " << setw(8) << result->delta_x_fix[1]<< " delta_Z_fix= " << setw(8) << result->delta_x_fix[2]<<
			"  dE(m)= " << setw(8) << result->dENU_fix[0] << "  dN(m)= " << setw(8) << result->dENU_fix[1] << "  dU(m)= " << setw(8) << result->dENU_fix[2] <<
			" delta_E_Fix= " << setw(8) << result->delta_e_fix[0]<< " delta_N_fix= " << setw(8) << result->delta_e_fix[1]<< " delta_U_fix= " << setw(8) << result->delta_e_fix[2]<<
			"  dX_Std(m)= " << setw(8) << result->dx_FixStd[0] << "  dY_Std(m)= " << setw(8) << result->dx_FixStd[1] << "  dZ_Std(m)= " << setw(8) << result->dx_FixStd[2] <<
			"  dE_Std(m)= " << setw(8) << result->dENU_FixStd[0] << "  dN_Std(m)= " << setw(8) << result->dENU_FixStd[1] << "  dU_Std(m)= " << setw(8) << result->dENU_FixStd[2] <<
			"  dX_Flo(m)= " << setw(8) << result->dX_flo[0] << "  dY_Flo(m)= " << setw(8) << result->dX_flo[1] << "  dZ_Flo(m)= " << setw(8) << result->dX_flo[2] <<
			" delta_X_flo= " << setw(8) << result->delta_x_flo[0]<< " delta_Y_flo= " << setw(8) << result->delta_x_flo[1]<< " delta_Z_flo= " << setw(8) << result->delta_x_flo[2]<<
			"  dE_Flo(m)= " << setw(8) << result->dENU_flo[0] << "  dN_Flo(m)= " << setw(8) << result->dENU_flo[1] << "  dU_Flo(m)= " << setw(8) << result->dENU_flo[2] <<
			" delta_E_flo= " << setw(8) << result->delta_e_flo[0]<< " delta_N_flo= " << setw(8) << result->delta_e_flo[1]<< " delta_U_flo= " << setw(8) << result->delta_e_flo[2]<<
			"  dX_Flo_Std(m)= " << setw(8) << result->mxyz[0] << "  dY_Flo_Std(m)= " << setw(8) << result->mxyz[1] << "  dZ_Flo_Std(m)= " << setw(8) << result->mxyz[2] <<
			"  dE_Flo_Std(m)= " << setw(8) << result->mENU[0] << "  dN_Flo_Std(m)= " << setw(8) << result->mENU[1] << "  dU_Flo_Std(m)= " << setw(8) << result->mENU[2] <<
			"  Sigma(m)= " << setw(8) << result->sigma << "  RDOP(m)= " << setw(8) << result->RDOP <<" RMS= " << setw(8) <<result->RMS << setprecision(3) << "  Ratio= " << setw(8) << result->Ratio <<
			"  GPSNum= " << result->GPSNum << "  BDSNum= " << result->BDSNum << endl;
	}
	else if (result->Ratio<=2.0)
	{
		cout << result->Time.Week << " " << fixed << setprecision(3) << result->Time.SecOfWeek << " Flo " << "GPSNum=" << result->GPSNum << " BDSNum=" << result->BDSNum << endl;
		outfile << result->Time.Week << " " << fixed << setprecision(3) << result->Time.SecOfWeek << " Flo " << " dX(m)= " <<
			fixed << setprecision(4) << setw(8) << result->dX_fix[0] << "  dY(m)= " << setw(8) << result->dX_fix[1] << "  dZ(m)= " << setw(8) << result->dX_fix[2] <<
			" delta_X_fix= " << setw(8) << result->delta_x_fix[0] << " delta_Y_fix= " << setw(8) << result->delta_x_fix[1] << " delta_Z_fix= " << setw(8) << result->delta_x_fix[2] <<
			"  dE(m)= " << setw(8) << result->dENU_fix[0] << "  dN(m)= " << setw(8) << result->dENU_fix[1] << "  dU(m)= " << setw(8) << result->dENU_fix[2] <<
			" delta_E_Fix= " << setw(8) << result->delta_e_fix[0] << " delta_N_fix= " << setw(8) << result->delta_e_fix[1] << " delta_U_fix= " << setw(8) << result->delta_e_fix[2] <<
			"  dX_Std(m)= " << setw(8) << result->dx_FixStd[0] << "  dY_Std(m)= " << setw(8) << result->dx_FixStd[1] << "  dZ_Std(m)= " << setw(8) << result->dx_FixStd[2] <<
			"  dE_Std(m)= " << setw(8) << result->dENU_FixStd[0] << "  dN_Std(m)= " << setw(8) << result->dENU_FixStd[1] << "  dU_Std(m)= " << setw(8) << result->dENU_FixStd[2] <<
			"  dX_Flo(m)= " << setw(8) << result->dX_flo[0] << "  dY_Flo(m)= " << setw(8) << result->dX_flo[1] << "  dZ_Flo(m)= " << setw(8) << result->dX_flo[2] <<
			" delta_X_flo= " << setw(8) << result->delta_x_flo[0] << " delta_Y_flo= " << setw(8) << result->delta_x_flo[1] << " delta_Z_flo= " << setw(8) << result->delta_x_flo[2] <<
			"  dE_Flo(m)= " << setw(8) << result->dENU_flo[0] << "  dN_Flo(m)= " << setw(8) << result->dENU_flo[1] << "  dU_Flo(m)= " << setw(8) << result->dENU_flo[2] <<
			" delta_E_flo= " << setw(8) << result->delta_e_flo[0] << " delta_N_flo= " << setw(8) << result->delta_e_flo[1] << " delta_U_flo= " << setw(8) << result->delta_e_flo[2] <<
			"  dX_Flo_Std(m)= " << setw(8) << result->mxyz[0] << "  dY_Flo_Std(m)= " << setw(8) << result->mxyz[1] << "  dZ_Flo_Std(m)= " << setw(8) << result->mxyz[2] <<
			"  dE_Flo_Std(m)= " << setw(8) << result->mENU[0] << "  dN_Flo_Std(m)= " << setw(8) << result->mENU[1] << "  dU_Flo_Std(m)= " << setw(8) << result->mENU[2] <<
			"  Sigma(m)= " << setw(8) << result->sigma << "  RDOP(m)= " << setw(8) << result->RDOP << " RMS= " << setw(8) << result->RMS << setprecision(3) << "  Ratio= " << setw(8) << result->Ratio <<
			"  GPSNum= " << result->GPSNum << "  BDSNum= " << result->BDSNum << endl;
	}
}