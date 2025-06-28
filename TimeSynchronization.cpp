#include "header.h"
#include "sockets.h"

/*
****************************************************************************
函数名：数据同步函数（事后）
参数：FBas    基站数据文件的文件流
      FRov    流动站数据文件的文件流
      rtkdata 同步数据存储处
返回值：数据同步结果的标志  1-数据同步成功  0-数据同步失败 -1-文件数据结束
函数功能：读取并解码两个站对应文件中的数据，并将其进行时间同步
****************************************************************************
*/
int TimeSyncPP(ifstream& FBas, ifstream& FRov, rtkData_t* rtkdata)
{
    // 将上一历元的两站数据存储好
    memcpy(&rtkdata->basData0, &rtkdata->basData, sizeof(rtkdata->basData0));
    memcpy(&rtkdata->rovData0, &rtkdata->rovData, sizeof(rtkdata->rovData0));

    //解码流动站数据
   static unsigned char Rbuf[MAXRAWLEN];  //存储解码数据包的数组
   static int ReadRLen = 0; //已读取的数据量
   static int Rjudge = 0;   //是否读取到流动站观测数据的标识（1---读取到，0---未读取到）
    while (!FRov.eof())
    {
        FRov.read((char*)(Rbuf + ReadRLen), MAXRAWLEN - ReadRLen);
        int Rlength = static_cast<int>(FRov.gcount()); //上一步read的字节数
        if (Rlength == 0) break;
        ReadRLen += Rlength;
        Rjudge = DecodeNovOem7Dat(Rbuf, ReadRLen, &rtkdata->rovData, rtkdata->GPSEph, rtkdata->BDSEph, &rtkdata->RovPres);
        if (Rjudge == 1) break; // 解码到了观测值，跳出循环
    }

   // if (rtkdata->rovData.Time.SecOfWeek == 11328)
    //{
    //    rtkdata->rovData.Time.SecOfWeek = 11328;
   // }
    // 求两站时间差值(单位为s)
    double dt = diffTime(&rtkdata->rovData.Time, &rtkdata->basData.Time);

    if (fabs(dt) <= 0.001) return 1; // 时间同步
    else if (dt < 0) return 0; // 流动站数据滞后，无法解算
    else // 基站数据滞后，需要读取更后面时刻的基站数据(循环执行直至时间同步为止)
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
                int Blength = static_cast<int>(FBas.gcount()); //上一步read的字节数
                if (Blength == 0) break;
                ReadBLen += Blength;
                Bjudge = DecodeNovOem7Dat(Bbuf, ReadBLen, &rtkdata->basData, rtkdata->GPSEph, rtkdata->BDSEph, &rtkdata->BasPres);
                if (Bjudge == 1) break; // 解码到了观测值，跳出循环
            }

            dt = diffTime(&rtkdata->rovData.Time, &rtkdata->basData.Time);

            if (fabs(dt) <= 0.001)  flag = false;// 时间完成同步
            else if (dt < 0) return 0; // 基站历元数据丢失，无法解算

        } while (flag);
        return 1;
    }
}

/*
****************************************************************************
函数名：数据同步函数（实时）
参数：BasSock    基站数据的网口
      RovSock    流动站数据的网口
      rtkdata    同步数据存储处（两个站的观测值、星历、参考定位结果）
      BasData    存储基准站原始二进制数据的文件流
      RovData    存储流动站原始二进制数据的文件流
返回值：数据同步结果的标志  1-数据同步成功  0-数据同步失败
函数功能：读取并解码两个站对应网口中的数据，并将其进行时间同步
****************************************************************************
*/
int TimeSyncRT(SOCKET& BasSock, SOCKET& RovSock, rtkData_t* rtkdata,vector<epoch_t>& rovOBS, vector<epoch_t>& basOBS,ofstream& BasData,ofstream& RovData)
{
    // 将上一历元的两站数据存储好
    memcpy(&rtkdata->basData0, &rtkdata->basData, sizeof(rtkdata->basData0));
    memcpy(&rtkdata->rovData0, &rtkdata->rovData, sizeof(rtkdata->rovData0));

    //Sleep(800);

    int irove, ibase;

    // 解码流动站数据
    static unsigned char RBuf[MAXRAWLEN];
    unsigned Rbuf[MAXBUFF];
    //static unsigned char RBuf[MAXRAWLEN * 10];
    //unsigned char Rbuf[MAXRAWLEN * 5];
    static int RlenB = 0; // RBuf的长度（字节数）
    static unsigned char BBuf[MAXRAWLEN];
    unsigned Bbuf[MAXBUFF];
    //static unsigned char BBuf[MAXRAWLEN*10];
    //unsigned char Bbuf[MAXRAWLEN*5];
    static int BlenB = 0;
    int BlenT = 0;
    int RlenT = 0; // 中间存储地的长度
    double dt;
    while (RovSock)
    {
        if ((RlenT = recv(RovSock, (char*)(RBuf + RlenB), MAXRAWLEN - RlenB, 0)) > 0)
        //if((RlenT = recv(RovSock, (char*)(Rbuf), MAXBUFF, 0))>0)
        {
            RovData.write((char*)(RBuf + RlenB), RlenT);//保存原始数据
            //RovData.write((char*)(Rbuf), RlenT);//保存原始数据
            //if ((RlenB + RlenT) > MAXRAWLEN)  
                //RlenB = 0; // 若本次读入数据加入大BUFF后会超过其最大长度，则先将原先已有的数据舍弃
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
                            BasData.write((char*)(BBuf + BlenB), BlenT);//保存原始数据
                            //BasData.write((char*)(Bbuf), BlenT);//保存原始数据
                            //if ((BlenB + BlenT) > MAXRAWLEN) 
                                //BlenB = 0; // 若本次读入数据加入大BUFF后会超过其最大长度，则先将原先已有的数据舍弃
                            //memcpy(BBuf + BlenB, Bbuf, BlenT);
                            BlenB = BlenT + BlenB;
                            if(DecodeNovOem7Dat(BBuf, BlenB, &rtkdata->basData, rtkdata->GPSEph, rtkdata->BDSEph, &rtkdata->BasPres)==1)
                            {
                                dt = diffTime(&rtkdata->rovData.Time, &rtkdata->basData.Time);
                                if (fabs(dt) < 0.2) return 1;
                                else if (dt < 0)
                                {
                                    cout << rtkdata->rovData.Time.Week << " " << rtkdata->rovData.Time.SecOfWeek << " 同步失败" << endl;
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
        RovData.write((char*)(Rbuf), RlenT);//保存原始数据
        if ((RlenB + RlenT) > MAXRAWLEN * 10)  RlenB = 0; // 若本次读入数据加入大BUFF后会超过其最大长度，则先将原先已有的数据舍弃
        memcpy(RBuf + RlenB, Rbuf, RlenT);
        RlenB = RlenT + RlenB;
        if (rovOBS.size() >= 5)
        {
            //cout << rovOBS[0].Time.Week << " " << rovOBS[0].Time.SecOfWeek << " 同步失败" << endl;
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
        BasData.write((char*)(Bbuf), BlenT);//保存原始数据
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