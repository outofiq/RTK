#include "header.h"

/*
*************************************************************************************
��������ѡȡ�ο���
������RovData   ָ������վ�۲����ݵ�ָ��
	  SDObs     ָ�򵥲����ݵ�ָ��
	  DDObs     ָ��˫�����ݵ�ָ��
�������ܣ�ѡȡ��׼�ǣ�˫ϵͳ��һ����������PRN�������������˫��������
ע�⣺�����RovData��Ҫ�ڱ���������ǰ����һ��SPP���õ�����λ�á��߶Ƚǵ��������
*************************************************************************************
*/
void SelectRefSat(const epoch_t* RovData, SDEpochOBS_t* SDObs, DDOBS_t* DDObs)
{
	short sys = -1; // 0---GPS,1---BDS
	double MaxEle[2] = { 0.0,0.0 }; //MaxEle[0]---GPS���߶Ƚǣ�MaxEle[1]---BDS���߶Ƚ�
	for (int i = 0; i < SDObs->SatNum; i++)
	{
		//ȷ������ϵͳ����
		if (SDObs->SdSatObs[i].System == sys_t::GPS) sys = 0;
		else if (SDObs->SdSatObs[i].System == sys_t::BDS) sys = 1;
		else continue;

		//������վ�ĸ߶Ƚ�Ϊ׼
		if (RovData->SatPVT[SDObs->SdSatObs[i].nRov].Elevation > MaxEle[sys])
		{
			MaxEle[sys] = RovData->SatPVT[SDObs->SdSatObs[i].nRov].Elevation;
			DDObs->RefPrn[sys] = SDObs->SdSatObs[i].Prn;
			DDObs->RefPos[sys] = i;
		}
		else continue;
	}
	if (DDObs->RefPos[0] > -1) DDObs->RefGPS = true;
	if (DDObs->RefPos[1] > -1) DDObs->RefBDS = true;
}

/*
*************************************************************************************
������������˫��۲�ֵ
������
	  SDObs     ָ�򵥲����ݵ�ָ��
	  DDObs     ָ��˫�����ݵ�ָ��
�������ܣ���������۲�ֵ������ÿ��������ο���֮���˫��۲�ֵ
ע�⣺�����DDObs��Ҫ�ڱ���������ǰѡȡ�ο���
*************************************************************************************
*/
void DDEpochObs(SDEpochOBS_t* SDObs, DDOBS_t* DDObs)
{
	int GPSDDOBSNum = 0;
	int BDSDDObsNum = 0;
	for (int i = 0; i < SDObs->SatNum; i++)
	{
		if (i == DDObs->RefPos[0] || i == DDObs->RefPos[1]) continue; // �ο����������޷����
		else;

		if (SDObs->SdSatObs[i].System == sys_t::GPS && DDObs->RefGPS)
		{
			DDObs->GPSDDP[0][GPSDDOBSNum] = SDObs->SdSatObs[i].dP[0] - SDObs->SdSatObs[DDObs->RefPos[0]].dP[0];
			DDObs->GPSDDP[1][GPSDDOBSNum] = SDObs->SdSatObs[i].dP[1] - SDObs->SdSatObs[DDObs->RefPos[0]].dP[1];
			DDObs->GPSDDL[0][GPSDDOBSNum] = SDObs->SdSatObs[i].dL[0] - SDObs->SdSatObs[DDObs->RefPos[0]].dL[0];
			DDObs->GPSDDL[1][GPSDDOBSNum] = SDObs->SdSatObs[i].dL[1] - SDObs->SdSatObs[DDObs->RefPos[0]].dL[1];
			DDObs->otherGPS[GPSDDOBSNum] = i;
			GPSDDOBSNum++;
		}
		else if (SDObs->SdSatObs[i].System == sys_t::BDS && DDObs->RefBDS)
		{
			DDObs->BDSDDP[0][BDSDDObsNum] = SDObs->SdSatObs[i].dP[0] - SDObs->SdSatObs[DDObs->RefPos[1]].dP[0];
			DDObs->BDSDDP[1][BDSDDObsNum] = SDObs->SdSatObs[i].dP[1] - SDObs->SdSatObs[DDObs->RefPos[1]].dP[1];
			DDObs->BDSDDL[0][BDSDDObsNum] = SDObs->SdSatObs[i].dL[0] - SDObs->SdSatObs[DDObs->RefPos[1]].dL[0];
			DDObs->BDSDDL[1][BDSDDObsNum] = SDObs->SdSatObs[i].dL[1] - SDObs->SdSatObs[DDObs->RefPos[1]].dL[1];
			DDObs->otherBDS[BDSDDObsNum] = i;
			BDSDDObsNum++;
		}
	}
	DDObs->DDSatNum[0] = GPSDDOBSNum;
	DDObs->DDSatNum[1] = BDSDDObsNum;
	DDObs->SatsNum = DDObs->DDSatNum[0] + DDObs->DDSatNum[1];
}