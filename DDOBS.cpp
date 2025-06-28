#include "header.h"

/*
*************************************************************************************
函数名：选取参考星
参数：RovData   指向流动站观测数据的指针
	  SDObs     指向单差数据的指针
	  DDObs     指向双差数据的指针
函数功能：选取基准星（双系统各一个），将其PRN号与索引存放在双差数据中
注意：传入的RovData需要在本函数外提前经过一次SPP来得到卫星位置、高度角等相关数据
*************************************************************************************
*/
void SelectRefSat(const epoch_t* RovData, SDEpochOBS_t* SDObs, DDOBS_t* DDObs)
{
	short sys = -1; // 0---GPS,1---BDS
	double MaxEle[2] = { 0.0,0.0 }; //MaxEle[0]---GPS最大高度角，MaxEle[1]---BDS最大高度角
	for (int i = 0; i < SDObs->SatNum; i++)
	{
		//确定卫星系统类型
		if (SDObs->SdSatObs[i].System == sys_t::GPS) sys = 0;
		else if (SDObs->SdSatObs[i].System == sys_t::BDS) sys = 1;
		else continue;

		//以流动站的高度角为准
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
函数名：计算双差观测值
参数：
	  SDObs     指向单差数据的指针
	  DDObs     指向双差数据的指针
函数功能：遍历单差观测值，计算每颗卫星与参考星之间的双差观测值
注意：传入的DDObs需要在本函数外提前选取参考星
*************************************************************************************
*/
void DDEpochObs(SDEpochOBS_t* SDObs, DDOBS_t* DDObs)
{
	int GPSDDOBSNum = 0;
	int BDSDDObsNum = 0;
	for (int i = 0; i < SDObs->SatNum; i++)
	{
		if (i == DDObs->RefPos[0] || i == DDObs->RefPos[1]) continue; // 参考卫星自身无法求差
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