#include "header.h"

/*
****************************************************************************
函数名：粗差探测函数
参数：	obs 观测数据
函数功能：对obs粗差探测并在SatObs里标记Valid,计算每颗卫星的（P）IF组合观测值
****************************************************************************
*/
void DetectError(epoch_t* obs)
{
	// 为了不破坏上一历元的组合观测值数组而建立的缓冲区
	MWGF_t comobs[MAXCHANNUM];// 存放当前历元所算组合值

	for (int i = 0; i < obs->SatNum; i++)
	{
		// 检查每一颗卫星的双频伪距和相位数据是否完整
		if (fabs(obs->SatObs[i].p[0]) < 1e-8 || fabs(obs->SatObs[i].p[1]) < 1e-8 || fabs(obs->SatObs[i].l[0]) < 1e-8 || fabs(obs->SatObs[i].l[1]) < 1e-8)
		{
			obs->SatObs[i].Valid = false;
			continue; //本颗卫星数据残缺，不计算组合观测值且观测值有效性设置为false
		}
		
		// 计算当前历元的该卫星的GF和MW组合值
		comobs[i].Sys = obs->SatObs[i].Sys;
		comobs[i].Prn = obs->SatObs[i].Prn;
		comobs[i].GF = obs->SatObs[i].l[0] - obs->SatObs[i].l[1];
		if (comobs[i].Sys == sys_t::GPS)
		{
			comobs[i].MW = (FG1_GPS * obs->SatObs[i].l[0] - FG2_GPS * obs->SatObs[i].l[1]) / (FG1_GPS - FG2_GPS) -
				(FG1_GPS * obs->SatObs[i].p[0] + FG2_GPS * obs->SatObs[i].p[1]) / (FG1_GPS + FG2_GPS);
			comobs[i].PIF = (FG1_GPS * FG1_GPS * obs->SatObs[i].p[0] - FG2_GPS * FG2_GPS * obs->SatObs[i].p[1]) / (FG1_GPS * FG1_GPS - FG2_GPS * FG2_GPS);
		}
		else if (comobs[i].Sys == sys_t::BDS)
		{
			comobs[i].MW = (FG1_BDS * obs->SatObs[i].l[0] - FG3_BDS * obs->SatObs[i].l[1]) / (FG1_BDS - FG3_BDS) -
				(FG1_BDS * obs->SatObs[i].p[0] + FG3_BDS * obs->SatObs[i].p[1]) / (FG1_BDS + FG3_BDS);
			comobs[i].PIF = (FG1_BDS * FG1_BDS * obs->SatObs[i].p[0] - FG3_BDS * FG3_BDS * obs->SatObs[i].p[1]) / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS);
		}

		comobs[i].n = 1;

		// 从上个历元的MWGF数据中查找该卫星的GF和MW组合值
		for (int j = 0; j < MAXCHANNUM; j++)
		{
			if (obs->ComObs[j].Sys == comobs[i].Sys && obs->ComObs[j].Prn == comobs[i].Prn)
			{
				// 找到了则比较二者是否在限差内
				double dGF = fabs(comobs[i].GF - obs->ComObs[j].GF);
				double dMW = fabs(comobs[i].MW - obs->ComObs[j].MW);

				// 没有超限，则标记true并计算MW的平滑值
				if (dGF < 0.05 && dMW < 3)
				{
					obs->SatObs[i].Valid = true;
					comobs[i].MW = (obs->ComObs[j].MW * obs->ComObs[j].n + comobs[i].MW) / (obs->ComObs[j].n + 1.0);
					comobs[i].n = obs->ComObs[j].n + 1;
				}
				else break;
				// 超限，则标记为粗差（而初始化本就是false，故不需处理）
			}
			// 如果在上个历元中没有找到，则不知其是否可用，默认为初始化时的false
			else continue;
		}
	}
	// 将缓冲区组合观测值的拷贝到obs里
	memcpy(obs->ComObs, comobs, sizeof(comobs));
}

/*
****************************************************
函数名：对流层改正函数
参数：hgt    测站高度m
	  elev   卫星相对于测站的高度角(rad)
返回值：对流层改正值
函数功能：根据测站高度和卫星高度角来算的对流层改正
****************************************************
*/
double hopfield(double hgt, double elev)
{
	double t0, p0, e0, h0;
	double t, p, e;
	double dend, elev2, denw, hw, hd, rkd, rkw;
	double trop;

	if (fabs(hgt) > 30000.0)   return 0.0;

	t0 = 20 + 273.16;
	p0 = 1013.0;
	e0 = 0.5 * exp(-37.2465 + 0.213166 * t0 - 0.000256908 * t0 * t0);
	h0 = 0;
	hw = 11000.0;
	t = t0 - 0.0068 * (hgt - h0);
	p = p0 * pow(1.0 - 0.0068 / t0 * (hgt - h0), 5);
	e = e0 * pow((1 - 0.0068 / t0 * (hgt - h0)), 2.0) * pow((1.0 - (hgt - h0) / hw), 4.0);
	elev2 = elev * elev * R2D * R2D;
	dend = sqrt(elev2 + 6.25) * D2R;
	denw = sqrt(elev2 + 2.25) * D2R;

	hd = 148.72 * t0 - 488.3552;
	rkd = 1.552e-5 * p / t * (hd - hgt);
	rkw = 7.46512e-2 * (e / t / t) * (hw - hgt);
	trop = (rkd / sin(dend)) + (rkw / sin(denw));
	return trop;
}