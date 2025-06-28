#include "header.h"
#include "sockets.h"

/*
****************************************************************************
������������ͬ���������º�
������FBas    ��վ�����ļ����ļ���
      FRov    ����վ�����ļ����ļ���
      rtkdata ͬ�����ݴ洢��
����ֵ������ͬ������ı�־  1-����ͬ���ɹ�  0-����ͬ��ʧ�� -1-�ļ����ݽ���
�������ܣ���ȡ����������վ��Ӧ�ļ��е����ݣ����������ʱ��ͬ��
****************************************************************************
*/
int TimeSyncPP(ifstream& FBas, ifstream& FRov, rtkData_t* rtkdata)
{
    // ����һ��Ԫ����վ���ݴ洢��
    memcpy(&rtkdata->basData0, &rtkdata->basData, sizeof(rtkdata->basData0));
    memcpy(&rtkdata->rovData0, &rtkdata->rovData, sizeof(rtkdata->rovData0));

    //��������վ����
   static unsigned char Rbuf[MAXRAWLEN];  //�洢�������ݰ�������
   static int ReadRLen = 0; //�Ѷ�ȡ��������
   static int Rjudge = 0;   //�Ƿ��ȡ������վ�۲����ݵı�ʶ��1---��ȡ����0---δ��ȡ����
    while (!FRov.eof())
    {
        FRov.read((char*)(Rbuf + ReadRLen), MAXRAWLEN - ReadRLen);
        int Rlength = static_cast<int>(FRov.gcount()); //��һ��read���ֽ���
        if (Rlength == 0) break;
        ReadRLen += Rlength;
        Rjudge = DecodeNovOem7Dat(Rbuf, ReadRLen, &rtkdata->rovData, rtkdata->GPSEph, rtkdata->BDSEph, &rtkdata->RovPres);
        if (Rjudge == 1) break; // ���뵽�˹۲�ֵ������ѭ��
    }

   // if (rtkdata->rovData.Time.SecOfWeek == 11328)
    //{
    //    rtkdata->rovData.Time.SecOfWeek = 11328;
   // }
    // ����վʱ���ֵ(��λΪs)
    double dt = diffTime(&rtkdata->rovData.Time, &rtkdata->basData.Time);

    if (fabs(dt) <= 0.001) return 1; // ʱ��ͬ��
    else if (dt < 0) return 0; // ����վ�����ͺ��޷�����
    else // ��վ�����ͺ���Ҫ��ȡ������ʱ�̵Ļ�վ����(ѭ��ִ��ֱ��ʱ��ͬ��Ϊֹ)
    {
        static unsigned char Bbuf[MAXRAWLEN];
        static int ReadBLen = 0;
        static int Bjudge = 0;
        bool flag = true;
        do
        {
            while (!FBas.eof())
            {
                FBas.read((char*)(Bbuf + ReadBLen), MAXRAWLEN - ReadBLen);
                int Blength = static_cast<int>(FBas.gcount()); //��һ��read���ֽ���
                if (Blength == 0) break;
                ReadBLen += Blength;
                Bjudge = DecodeNovOem7Dat(Bbuf, ReadBLen, &rtkdata->basData, rtkdata->GPSEph, rtkdata->BDSEph, &rtkdata->BasPres);
                if (Bjudge == 1) break; // ���뵽�˹۲�ֵ������ѭ��
            }

            dt = diffTime(&rtkdata->rovData.Time, &rtkdata->basData.Time);

            if (fabs(dt) <= 0.001)  flag = false;// ʱ�����ͬ��
            else if (dt < 0) return 0; // ��վ��Ԫ���ݶ�ʧ���޷�����

        } while (flag);
        return 1;
    }
}

/*
****************************************************************************
������������ͬ��������ʵʱ��
������BasSock    ��վ���ݵ�����
      RovSock    ����վ���ݵ�����
      rtkdata    ͬ�����ݴ洢��������վ�Ĺ۲�ֵ���������ο���λ�����
      BasData    �洢��׼վԭʼ���������ݵ��ļ���
      RovData    �洢����վԭʼ���������ݵ��ļ���
����ֵ������ͬ������ı�־  1-����ͬ���ɹ�  0-����ͬ��ʧ��
�������ܣ���ȡ����������վ��Ӧ�����е����ݣ����������ʱ��ͬ��
****************************************************************************
*/
int TimeSyncRT(SOCKET& BasSock, SOCKET& RovSock, rtkData_t* rtkdata,vector<epoch_t>& rovOBS, vector<epoch_t>& basOBS,ofstream& BasData,ofstream& RovData)
{
    // ����һ��Ԫ����վ���ݴ洢��
    memcpy(&rtkdata->basData0, &rtkdata->basData, sizeof(rtkdata->basData0));
    memcpy(&rtkdata->rovData0, &rtkdata->rovData, sizeof(rtkdata->rovData0));

    //Sleep(800);

    int irove, ibase;

    // ��������վ����
    static unsigned char RBuf[MAXRAWLEN];
    unsigned Rbuf[MAXBUFF];
    //static unsigned char RBuf[MAXRAWLEN * 10];
    //unsigned char Rbuf[MAXRAWLEN * 5];
    static int RlenB = 0; // RBuf�ĳ��ȣ��ֽ�����
    static unsigned char BBuf[MAXRAWLEN];
    unsigned Bbuf[MAXBUFF];
    //static unsigned char BBuf[MAXRAWLEN*10];
    //unsigned char Bbuf[MAXRAWLEN*5];
    static int BlenB = 0;
    int BlenT = 0;
    int RlenT = 0; // �м�洢�صĳ���
    double dt;
    while (RovSock)
    {
        if ((RlenT = recv(RovSock, (char*)(RBuf + RlenB), MAXRAWLEN - RlenB, 0)) > 0)
        //if((RlenT = recv(RovSock, (char*)(Rbuf), MAXBUFF, 0))>0)
        {
            RovData.write((char*)(RBuf + RlenB), RlenT);//����ԭʼ����
            //RovData.write((char*)(Rbuf), RlenT);//����ԭʼ����
            //if ((RlenB + RlenT) > MAXRAWLEN)  
                //RlenB = 0; // �����ζ������ݼ����BUFF��ᳬ������󳤶ȣ����Ƚ�ԭ�����е���������
            //memcpy(RBuf + RlenB, Rbuf, RlenT);
            RlenB = RlenT + RlenB;
            if (DecodeNovOem7Dat(RBuf, RlenB, &rtkdata->rovData, rtkdata->GPSEph, rtkdata->BDSEph, &rtkdata->RovPres) == 1)
            {
                dt= diffTime(&rtkdata->rovData.Time, &rtkdata->basData.Time);
                if (fabs(dt) < 0.2) return 1;
                else if (dt < 0) continue;
                else
                {
                    while (BasSock)
                    {
                        if ((BlenT = recv(BasSock, (char*)(BBuf + BlenB), MAXRAWLEN - BlenB, 0)) > 0)
                        //if ((BlenT = recv(BasSock, (char*)(Bbuf), MAXBUFF, 0)) > 0)
                        {
                            BasData.write((char*)(BBuf + BlenB), BlenT);//����ԭʼ����
                            //BasData.write((char*)(Bbuf), BlenT);//����ԭʼ����
                            //if ((BlenB + BlenT) > MAXRAWLEN) 
                                //BlenB = 0; // �����ζ������ݼ����BUFF��ᳬ������󳤶ȣ����Ƚ�ԭ�����е���������
                            //memcpy(BBuf + BlenB, Bbuf, BlenT);
                            BlenB = BlenT + BlenB;
                            if(DecodeNovOem7Dat(BBuf, BlenB, &rtkdata->basData, rtkdata->GPSEph, rtkdata->BDSEph, &rtkdata->BasPres)==1)
                            {
                                dt = diffTime(&rtkdata->rovData.Time, &rtkdata->basData.Time);
                                if (fabs(dt) < 0.2) return 1;
                                else if (dt < 0)
                                {
                                    cout << rtkdata->rovData.Time.Week << " " << rtkdata->rovData.Time.SecOfWeek << " ͬ��ʧ��" << endl;
                                    return 0;
                                }
                                else continue;
                            }
                        }
                    }
                    
                }
            }
        }
    }

    /*
    if ((RlenT = recv(RovSock, (char*)Rbuf, 5*MAXRAWLEN, 0)) > 0)
    {
        RovData.write((char*)(Rbuf), RlenT);//����ԭʼ����
        if ((RlenB + RlenT) > MAXRAWLEN * 10)  RlenB = 0; // �����ζ������ݼ����BUFF��ᳬ������󳤶ȣ����Ƚ�ԭ�����е���������
        memcpy(RBuf + RlenB, Rbuf, RlenT);
        RlenB = RlenT + RlenB;
        if (rovOBS.size() >= 5)
        {
            //cout << rovOBS[0].Time.Week << " " << rovOBS[0].Time.SecOfWeek << " ͬ��ʧ��" << endl;
            rovOBS.erase(rovOBS.begin());
        }
        if (DecodeNovOem7Dat(RBuf, RlenB, &rtkdata->rovData, rtkdata->GPSEph, rtkdata->BDSEph, &rtkdata->RovPres) == 1)
        {
            rovOBS.push_back(rtkdata->rovData);
            memset(&rtkdata->rovData, 0, sizeof(epoch_t));
        }
    }

    if ((BlenT = recv(BasSock, (char*)Bbuf, 5*MAXRAWLEN, 0)) > 0)
    {
        BasData.write((char*)(Bbuf), BlenT);//����ԭʼ����
        if ((BlenB + BlenT) > MAXRAWLEN * 10)  BlenB = 0;
        memcpy(BBuf + BlenB, Bbuf, BlenT);
        BlenB = BlenT + BlenB;
        if (basOBS.size() >= 5) basOBS.erase(basOBS.begin());
        if (DecodeNovOem7Dat(BBuf, BlenB, &rtkdata->basData, rtkdata->GPSEph, rtkdata->BDSEph, &rtkdata->BasPres)==1)
        {
            basOBS.push_back(rtkdata->basData);
            memset(&rtkdata->basData, 0, sizeof(epoch_t));
        }
    }

    timeSyn(rovOBS, basOBS, irove, ibase);
    if (irove == -1 || ibase == -1) return 0;
    rtkdata->rovData = rovOBS[irove];
    //rovOBS.erase(rovOBS.begin() + irove);
    rtkdata->basData = basOBS[ibase];
    return 1;
    */
}

void timeSyn(vector<epoch_t>rOBS, vector<epoch_t>bOBS, int& irove, int& ibase)
{
    char brove = 1, bbase = 1;

    irove = rOBS.size();// 5
    ibase = bOBS.size();// 5 

    while (1)
    {
        if (brove)irove--;
        if (bbase)ibase--;
        if (irove == -1 || ibase == -1)break;
        if (rOBS[irove].Time.SecOfWeek < 10 || bOBS[ibase].Time.SecOfWeek < 10)continue;

        if (fabs(diffTime(&rOBS[irove].Time, &bOBS[ibase].Time)) < 0.2)
            break;
        else if (rOBS[irove].Time.SecOfWeek < bOBS[ibase].Time.SecOfWeek)
        {
            brove = 0; bbase = 1; continue;
        }
        else
        {
            brove = 1; bbase = 0; continue;
        }
    }
}