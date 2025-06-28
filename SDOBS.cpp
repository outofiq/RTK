#include "header.h"

/*
**************************************************************
��������վ�䵥���
������RovData  ָ������վ�۲����ݵ�ָ��
      BasData  ָ���վ�۲����ݵ�ָ��
      SDObs    �������洢��
�������ܣ�������վ���վ�Ĺ۲�����������
**************************************************************
*/
void SDEpochObs(epoch_t* RovData, epoch_t* BasData, SDEpochOBS_t* SDObs)
{
    memset(SDObs->SdSatObs, 0.0, sizeof(SDObs->SdSatObs));

    SDObs->Time = RovData->Time;

    int satnum = 0; // ����۲����ݵĿ���������
    int GPSnum = 0; // ����۲����ݵĿ���GPS������
    int BDSnum = 0; // ����۲����ݵĿ���BDS������

    //�����С��30dBHZ���߶Ƚ�С��10���������Ϊ������
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

    // ����������վ�Ĺ۲�����
    for (int i = 0; i < RovData->SatNum; i++)
    {
        // �����ջ����ڲ�����ָ�겻�ϸ���������������������������������ݣ����������
        if (RovData->SatObs[i].Valid == false || RovData->SatPVT[i].Valid == false) continue;
        else;
        // �ڲ������վ�Ĺ۲�����
        for (int j = 0; j < BasData->SatNum; j++)
        {
            // �����ջ����ڲ�����ָ�겻�ϸ���������������������������������ݣ����������
            if (BasData->SatObs[j].Valid == false || BasData->SatPVT[j].Valid == false) continue;
            // ��û�����������ƥ��
            else if (BasData->SatObs[j].Sys == RovData->SatObs[i].Sys && BasData->SatObs[j].Prn == RovData->SatObs[i].Prn)
            {
                SDObs->SdSatObs[satnum].System = RovData->SatObs[i].Sys;
                SDObs->SdSatObs[satnum].Prn = RovData->SatObs[i].Prn;
                SDObs->SdSatObs[satnum].nRov = i;
                SDObs->SdSatObs[satnum].nBas = j;

                //����˫Ƶ�ĵ���۲�ֵ
                for (int k = 0; k < 2; k++)
                {
                    // ֻ������վ��α�඼��ֵ�����䵥�����˵�����ݶ�ʧ->�����㴦��
                    if (fabs(RovData->SatObs[i].p[k]) > 1e-8 && fabs(BasData->SatObs[j].p[k]) > 1e-8)
                    {
                        SDObs->SdSatObs[satnum].dP[k] = RovData->SatObs[i].p[k] - BasData->SatObs[j].p[k];
                    }
                    else SDObs->SdSatObs[satnum].dP[k] = 0;
                    // ֻ������վ����λ����ֵ�����䵥�����˵�����ݶ�ʧ->�����㴦��
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