#include "header.h"

/*
****************************************************************************
���������ֲ�̽�⺯��
������	obs �۲�����
�������ܣ���obs�ֲ�̽�Ⲣ��SatObs����Valid,����ÿ�����ǵģ�P��IF��Ϲ۲�ֵ
****************************************************************************
*/
void DetectError(epoch_t* obs)
{
	// Ϊ�˲��ƻ���һ��Ԫ����Ϲ۲�ֵ����������Ļ�����
	MWGF_t comobs[MAXCHANNUM];// ��ŵ�ǰ��Ԫ�������ֵ

	for (int i = 0; i < obs->SatNum; i++)
	{
		// ���ÿһ�����ǵ�˫Ƶα�����λ�����Ƿ�����
		if (fabs(obs->SatObs[i].p[0]) < 1e-8 || fabs(obs->SatObs[i].p[1]) < 1e-8 || fabs(obs->SatObs[i].l[0]) < 1e-8 || fabs(obs->SatObs[i].l[1]) < 1e-8)
		{
			obs->SatObs[i].Valid = false;
			continue; //�����������ݲ�ȱ����������Ϲ۲�ֵ�ҹ۲�ֵ��Ч������Ϊfalse
		}
		
		// ���㵱ǰ��Ԫ�ĸ����ǵ�GF��MW���ֵ
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

		// ���ϸ���Ԫ��MWGF�����в��Ҹ����ǵ�GF��MW���ֵ
		for (int j = 0; j < MAXCHANNUM; j++)
		{
			if (obs->ComObs[j].Sys == comobs[i].Sys && obs->ComObs[j].Prn == comobs[i].Prn)
			{
				// �ҵ�����Ƚ϶����Ƿ����޲���
				double dGF = fabs(comobs[i].GF - obs->ComObs[j].GF);
				double dMW = fabs(comobs[i].MW - obs->ComObs[j].MW);

				// û�г��ޣ�����true������MW��ƽ��ֵ
				if (dGF < 0.05 && dMW < 3)
				{
					obs->SatObs[i].Valid = true;
					comobs[i].MW = (obs->ComObs[j].MW * obs->ComObs[j].n + comobs[i].MW) / (obs->ComObs[j].n + 1.0);
					comobs[i].n = obs->ComObs[j].n + 1;
				}
				else break;
				// ���ޣ�����Ϊ�ֲ����ʼ��������false���ʲ��账��
			}
			// ������ϸ���Ԫ��û���ҵ�����֪���Ƿ���ã�Ĭ��Ϊ��ʼ��ʱ��false
			else continue;
		}
	}
	// ����������Ϲ۲�ֵ�Ŀ�����obs��
	memcpy(obs->ComObs, comobs, sizeof(comobs));
}

/*
****************************************************
���������������������
������hgt    ��վ�߶�m
	  elev   ��������ڲ�վ�ĸ߶Ƚ�(rad)
����ֵ�����������ֵ
�������ܣ����ݲ�վ�߶Ⱥ����Ǹ߶Ƚ�����Ķ��������
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