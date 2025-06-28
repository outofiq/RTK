#include "header.h"

/*
**************************************************************
函数名：站间单差函数
参数：RovData  指向流动站观测数据的指针
      BasData  指向基站观测数据的指针
      SDObs    单差结果存储处
函数功能：将流动站与基站的观测数据求差并储存
**************************************************************
*/
void SDEpochObs(epoch_t* RovData, epoch_t* BasData, SDEpochOBS_t* SDObs)
{
    memset(SDObs->SdSatObs, 0.0, sizeof(SDObs->SdSatObs));

    SDObs->Time = RovData->Time;

    int satnum = 0; // 单差观测数据的可用卫星数
    int GPSnum = 0; // 单差观测数据的可用GPS卫星数
    int BDSnum = 0; // 单差观测数据的可用BDS卫星数

    //载噪比小于30dBHZ，高度角小于10°的卫星置为不可用
    for (int i = 0; i < RovData->SatNum; i++)
    {
        if (RovData->SatObs[i].CN0[0] <= 30 || RovData->SatObs[i].CN0[1] <= 30 || RovData->SatPVT[i].Elevation <= (10 / 180) * PI)
        {
            RovData->SatObs[i].Valid = false;
            RovData->SatPVT[i].Valid = false;
        }
    }
    for (int i = 0; i < BasData->SatNum; i++)
    {
        if (BasData->SatObs[i].CN0[0] <= 30 || BasData->SatObs[i].CN0[1] <= 30 || BasData->SatPVT[i].Elevation <= (10 / 180) * PI)
        {
            BasData->SatObs[i].Valid = false;
            BasData->SatPVT[i].Valid = false;
        }
    }

    // 外层遍历流动站的观测数据
    for (int i = 0; i < RovData->SatNum; i++)
    {
        // 若接收机的内部质量指标不合格或卫星星历不可用则跳过本颗卫星数据，不做单差处理
        if (RovData->SatObs[i].Valid == false || RovData->SatPVT[i].Valid == false) continue;
        else;
        // 内层遍历基站的观测数据
        for (int j = 0; j < BasData->SatNum; j++)
        {
            // 若接收机的内部质量指标不合格或卫星星历不可用则跳过本颗卫星数据，不做单差处理
            if (BasData->SatObs[j].Valid == false || BasData->SatPVT[j].Valid == false) continue;
            // 若没有问题则进行匹配
            else if (BasData->SatObs[j].Sys == RovData->SatObs[i].Sys && BasData->SatObs[j].Prn == RovData->SatObs[i].Prn)
            {
                SDObs->SdSatObs[satnum].System = RovData->SatObs[i].Sys;
                SDObs->SdSatObs[satnum].Prn = RovData->SatObs[i].Prn;
                SDObs->SdSatObs[satnum].nRov = i;
                SDObs->SdSatObs[satnum].nBas = j;

                //计算双频的单差观测值
                for (int k = 0; k < 2; k++)
                {
                    // 只有两个站的伪距都有值才做其单差，否则说明数据丢失->做置零处理
                    if (fabs(RovData->SatObs[i].p[k]) > 1e-8 && fabs(BasData->SatObs[j].p[k]) > 1e-8)
                    {
                        SDObs->SdSatObs[satnum].dP[k] = RovData->SatObs[i].p[k] - BasData->SatObs[j].p[k];
                    }
                    else SDObs->SdSatObs[satnum].dP[k] = 0;
                    // 只有两个站的相位都有值才做其单差，否则说明数据丢失->做置零处理
                    if (fabs(RovData->SatObs[i].l[k]) > 1e-8 && fabs(BasData->SatObs[j].l[k]) > 1e-8)
                    {
                        SDObs->SdSatObs[satnum].dL[k] = RovData->SatObs[i].l[k] - BasData->SatObs[j].l[k];
                    }
                    else SDObs->SdSatObs[satnum].dL[k] = 0;
                }

                if (RovData->SatObs[i].Sys == sys_t::GPS) GPSnum++;
                else if (RovData->SatObs[i].Sys == sys_t::BDS) BDSnum++;
                else;

                satnum++;
                break;
            }
            else continue;
        }
    }
    SDObs->SatNum = satnum;
    SDObs->GPSNum = GPSnum;
    SDObs->BDSNum = BDSnum;
}