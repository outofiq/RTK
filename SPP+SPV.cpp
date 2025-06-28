#include "header.h"

/*
**************************************************************************
函数名：单点定位
参数：obs     当前历元的观测数据
      GPSEph  GPS星历
      BDSEph  BDS星历
      pos     用户定位结果
             （传入时为上一历元的定位结果，函数运行后为本次历元的定位结果）
返回值：单点定位是否完成 true-成功 false-失败
函数功能：进行SPP解算并将结果存入pos
**************************************************************************
*/
bool SPP(epoch_t* obs, gpseph_t* GPSEph, gpseph_t* BDSEph, epochPOS_t* pos,CFGINFO_t cfg)
{
    pos->Time = obs->Time;
    // 设定初始位置
    // 第一个历元的初始位置设置为0（构造函数已做）
    // 后续历元的初始位置设置为上一历元解算结果
    XYZ X0;
    X0.x = pos->Pos[0];
    X0.y = pos->Pos[1];
    X0.z = pos->Pos[2];
    double dt0[2] = { 0,0 };// 接收机钟差（0—GPS，1—BDS）
    // 迭代标志
    bool flag = true;
    // 迭代计数器
    int calcu_num = 0;
    // 迭代解算，直至收敛
    do 
    {
        // 计算信号发射时刻的卫星位置、钟差、地球自转改正和对流层延迟
        CalSatPos(obs, GPSEph, BDSEph, X0);
        // 对所有卫星的观测数据进行线性化
        // 以初始位置为参考，对观测方程线性化，得B和W矩阵，统计参与定位的各系统卫星数和所有卫星数
        int satnum[2] = { 0,0 };// 第一个元素为GPS可用卫星数，第二个元素为BDS可用卫星数
         // 总的可用卫星数计数
        int sum_satnum = 0;

        //只计算GPS系统
        if (cfg.systype == 0)
        {
            Matrix B(MAXGPSPRN, 4, 0);
            Matrix w(MAXGPSPRN, 1, 0);
            Matrix P(MAXGPSPRN, MAXGPSPRN, 0);
            for (int i = 0; i < obs->SatNum; i++)
            {
                // 观测数据不完整或有粗差、卫星位置计算失败，不参与定位计算
                if (obs->SatObs[i].Valid == false || obs->SatPVT[i].Valid == false) continue;
                if (obs->SatObs[i].Sys != sys_t::GPS) continue;
                double ρ = sqrt((obs->SatPVT[i].Pos[0] - X0.x) * (obs->SatPVT[i].Pos[0] - X0.x) +
                    (obs->SatPVT[i].Pos[1] - X0.y) * (obs->SatPVT[i].Pos[1] - X0.y) +
                    (obs->SatPVT[i].Pos[2] - X0.z) * (obs->SatPVT[i].Pos[2] - X0.z));
                B(satnum[0], 0) = (X0.x - obs->SatPVT[i].Pos[0]) / ρ;
                B(satnum[0], 1) = (X0.y - obs->SatPVT[i].Pos[1]) / ρ;
                B(satnum[0], 2) = (X0.z - obs->SatPVT[i].Pos[2]) / ρ;
                B(satnum[0], 3) = 1;
                w(satnum[0], 0) = obs->ComObs[i].PIF - (ρ + dt0[0] - CLIGHT * obs->SatPVT[i].ClkOft + obs->SatPVT[i].TropCorr);
                satnum[0] ++;
                P(satnum[0], satnum[0]) = 1;
            }
            int x_num = 4;// 参数个数
            // 可用卫星数不足以解算，则退出函数
            if (satnum[0] < x_num) return false;
            Matrix N = B.transpose() * P * B;
            Matrix W = B.transpose() * P * w;
            Matrix x(4, 1, 0);
            Matrix N_inv = N.inverse();
            x = N_inv * W;
            double x_norm = sqrt(x(0, 0) * x(0, 0) + x(1, 0) * x(1, 0) + x(2, 0) * x(2, 0) + x(3, 0) * x(3, 0));
            if (x_norm < 1e-5)  flag = false;
            X0.x = X0.x + x(0, 0);
            X0.y = X0.y + x(1, 0);
            X0.z = X0.z + x(2, 0);
            // 定位精度评价，计算PDOP
            pos->PDOP = sqrt(N_inv(0, 0) + N_inv(1, 1) + N_inv(2, 2));
            dt0[0] = x(3, 0);
            // 定位精度评价，计算验后单位权中误差
            Matrix v = B * x - w;
            Matrix vTv = v.transpose() * P * v;
            pos->SigmaPos = sqrt(vTv(0, 0) / (satnum[0] - x_num));
            calcu_num++;
            pos->SatNum = satnum[0];
            if (calcu_num > 15) flag = false;
        }

        //只计算BDS系统
        if (cfg.systype == 1)
        {
            Matrix B(MAXBDSPRN, 4, 0);
            Matrix w(MAXBDSPRN, 1, 0);
            Matrix P(MAXBDSPRN, MAXBDSPRN, 0);
            for (int i = 0; i < obs->SatNum; i++)
            {
                // 观测数据不完整或有粗差、卫星位置计算失败，不参与定位计算
                if (obs->SatObs[i].Valid == false || obs->SatPVT[i].Valid == false) continue;
                if (obs->SatObs[i].Sys != sys_t::BDS) continue;
                double ρ = sqrt((obs->SatPVT[i].Pos[0] - X0.x) * (obs->SatPVT[i].Pos[0] - X0.x) +
                    (obs->SatPVT[i].Pos[1] - X0.y) * (obs->SatPVT[i].Pos[1] - X0.y) +
                    (obs->SatPVT[i].Pos[2] - X0.z) * (obs->SatPVT[i].Pos[2] - X0.z));
                B(satnum[1], 0) = (X0.x - obs->SatPVT[i].Pos[0]) / ρ;
                B(satnum[1], 1) = (X0.y - obs->SatPVT[i].Pos[1]) / ρ;
                B(satnum[1], 2) = (X0.z - obs->SatPVT[i].Pos[2]) / ρ;
                B(satnum[1], 3) = 1;
                w(satnum[1], 0) = obs->ComObs[i].PIF - (ρ + dt0[1] - CLIGHT * obs->SatPVT[i].ClkOft + obs->SatPVT[i].TropCorr + CLIGHT * (FG1_BDS * FG1_BDS * obs->SatPVT[i].Tgd1) / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS));
                satnum[1] ++;
                P(satnum[1], satnum[1]) = 1;
            }
            int x_num = 4;// 参数个数
            // 可用卫星数不足以解算，则退出函数
            if (satnum[1] < x_num) return false;
            Matrix N = B.transpose() * P * B;
            Matrix W = B.transpose() * P * w;
            Matrix x(4, 1, 0);
            Matrix N_inv = N.inverse();
            x = N_inv * W;
            double x_norm = sqrt(x(0, 0) * x(0, 0) + x(1, 0) * x(1, 0) + x(2, 0) * x(2, 0) + x(3, 0) * x(3, 0));
            if (x_norm < 1e-5)  flag = false;
            X0.x = X0.x + x(0, 0);
            X0.y = X0.y + x(1, 0);
            X0.z = X0.z + x(2, 0);
            // 定位精度评价，计算PDOP
            pos->PDOP = sqrt(N_inv(0, 0) + N_inv(1, 1) + N_inv(2, 2));
            dt0[1] = x(3, 0);
            // 定位精度评价，计算验后单位权中误差
            Matrix v = B * x - w;
            Matrix vTv = v.transpose() * P * v;
            pos->SigmaPos = sqrt(vTv(0, 0) / (satnum[1] - x_num));
            calcu_num++;
            pos->SatNum = satnum[1];
            if (calcu_num > 15) flag = false;
        }
    
        //GPS+BDS双系统
        if (cfg.systype == 2)
        {
            Matrix B(MAXCHANNUM, 5, 0);
            Matrix w(MAXCHANNUM, 1, 0);
            // 权阵
            Matrix P(MAXCHANNUM, MAXCHANNUM, 0);

            for (int i = 0; i < obs->SatNum; i++)
            {
                sum_satnum = satnum[0] + satnum[1];
                // 观测数据不完整或有粗差、卫星位置计算失败，不参与定位计算
                if (obs->SatObs[i].Valid == false || obs->SatPVT[i].Valid == false) continue;
                double ρ = sqrt((obs->SatPVT[i].Pos[0] - X0.x) * (obs->SatPVT[i].Pos[0] - X0.x) +
                    (obs->SatPVT[i].Pos[1] - X0.y) * (obs->SatPVT[i].Pos[1] - X0.y) +
                    (obs->SatPVT[i].Pos[2] - X0.z) * (obs->SatPVT[i].Pos[2] - X0.z));

                B(sum_satnum, 0) = (X0.x - obs->SatPVT[i].Pos[0]) / ρ;
                B(sum_satnum, 1) = (X0.y - obs->SatPVT[i].Pos[1]) / ρ;
                B(sum_satnum, 2) = (X0.z - obs->SatPVT[i].Pos[2]) / ρ;
                if (obs->SatObs[i].Sys == sys_t::GPS)
                {
                    B(sum_satnum, 3) = 1;
                    B(sum_satnum, 4) = 0;
                    w(sum_satnum, 0) = obs->ComObs[i].PIF - (ρ + dt0[0] - CLIGHT * obs->SatPVT[i].ClkOft + obs->SatPVT[i].TropCorr);
                    //w(sum_satnum, 0) = obs->ComObs[i].PIF - (ρ + dt0[0] - CLIGHT * obs->SatPVT[i].ClkOft ); //不进行对流层改正
                    satnum[0]++;// GPS可用卫星数加一
                }
                else if (obs->SatObs[i].Sys == sys_t::BDS)
                {

                    B(sum_satnum, 3) = 0;
                    B(sum_satnum, 4) = 1;
                    // 注意BDS的IF组合的观测方程会多一个tgd硬件延迟
                    w(sum_satnum, 0) = obs->ComObs[i].PIF - (ρ + dt0[1] - CLIGHT * obs->SatPVT[i].ClkOft + obs->SatPVT[i].TropCorr + CLIGHT * (FG1_BDS * FG1_BDS * obs->SatPVT[i].Tgd1) / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS));
                    //w(sum_satnum, 0) = obs->ComObs[i].PIF - (ρ + dt0[1] - CLIGHT * obs->SatPVT[i].ClkOft + CLIGHT * (FG1_BDS * FG1_BDS * obs->SatPVT[i].Tgd1) / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS)); //不进行对流层改正
                    satnum[1]++;// BDS可用卫星数加一
                }
                P(sum_satnum, sum_satnum) = 1;// 等权法
            }

            sum_satnum = satnum[0] + satnum[1];
            // 如果卫星数不足，直接返回定位失败
            int x_num = 3;// 参数个数
            for (int i = 0; i < 2; i++)
            {
                if (satnum[i] != 0) x_num++;
            }
            // 可用卫星数不足以解算，则退出函数
            if (sum_satnum < x_num) return false;
            // 建立法方程N，W
            Matrix N = B.transpose() * P * B;
            Matrix W = B.transpose() * P * w;
            Matrix x(5, 1, 0);

            // GPS\BDS双系统均有数据
            if (satnum[0] != 0 && satnum[1] != 0)
            {
                // 最小二乘求解 
                Matrix N_inv = N.inverse();
                x = N_inv * W;
                double x_norm = sqrt(x(0, 0) * x(0, 0) + x(1, 0) * x(1, 0) + x(2, 0) * x(2, 0) + x(3, 0) * x(3, 0) + x(4, 0) * x(4, 0));//sign
                if (x_norm < 1e-5)  flag = false;
                //更新定位结果
                X0.x = X0.x + x(0, 0);
                X0.y = X0.y + x(1, 0);
                X0.z = X0.z + x(2, 0);
                // 定位精度评价，计算PDOP
                pos->PDOP = sqrt(N_inv(0, 0) + N_inv(1, 1) + N_inv(2, 2));
                //更新钟差
                dt0[0] = x(3, 0);
                dt0[1] = x(4, 0);
            }
            // 只有GPS或BDS单个系统的数据
            else
            {
                // 重构N，W
                if (satnum[0] == 0)
                {
                    // 如果GPS卫星数为0，则第4行和第4列均为0，可以删除，将矩阵缩减为4 * 4
                    N.deletRow(3);
                    N.deleteColumn(3);
                    // 将W矩阵删除第4行
                    W.deletRow(3);
                    //将B矩阵删除第4列
                    B.deleteColumn(3);
                }
                else if (satnum[1] == 0)
                {
                    // BDS卫星数为0，将第5行和第5列删除
                    N.deletRow(4);
                    N.deleteColumn(4);
                    // 将W矩阵删除第5行
                    W.deletRow(4);
                    //将B矩阵删除第5列
                    B.deleteColumn(4);
                }
                // 最小二乘求解
                Matrix N_inv = N.inverse();
                x = N_inv * W;

                double x_norm = sqrt(x(0, 0) * x(0, 0) + x(1, 0) * x(1, 0) + x(2, 0) * x(2, 0) + x(3, 0) * x(3, 0));
                if (x_norm < 1e-5)  flag = false;
                X0.x = X0.x + x(0, 0);
                X0.y = X0.y + x(1, 0);
                X0.z = X0.z + x(2, 0);
                // 定位精度评价，计算PDOP
                pos->PDOP = sqrt(N_inv(0, 0) + N_inv(1, 1) + N_inv(2, 2));

                //更新钟差
                if (satnum[0] == 0)
                {
                    dt0[0] = 0;
                    dt0[1] = x(3, 0);
                }
                else if (satnum[1] == 0)
                {
                    dt0[0] = x(3, 0);
                    dt0[1] = 0;
                }
            }

            // 定位精度评价，计算验后单位权中误差
            Matrix v = B * x - w;
            Matrix vTv = v.transpose() * P * v;
            pos->SigmaPos = sqrt(vTv(0, 0) / (sum_satnum - x_num));
            calcu_num++;
            pos->SatNum = sum_satnum;
            if (calcu_num > 15) flag = false;
        }
    } while (flag);
    
    pos->Pos[0] = X0.x;
    pos->Pos[1] = X0.y;
    pos->Pos[2] = X0.z;

    return true;
}

/*
***************************************
函数名：单点测速
参数：obs     当前历元的观测数据
      pos     用户定位结果
函数功能：进行SPV解算并将结果存入pos
注：只有当SPP成功才能进行SPV
***************************************
*/
void SPV(epoch_t* obs, epochPOS_t* pos)
{
    Matrix B(MAXCHANNUM, 4, 0);
    Matrix w(MAXCHANNUM, 1, 0);
    int valid_satnum = 0;
    for (int i = 0; i < obs->SatNum; i++)
    {
        // 观测数据不完整或有粗差、卫星位置计算失败，不参与定位计算
        if (obs->SatObs[i].Valid == false || obs->SatPVT[i].Valid == false) continue;
        //if (obs->SatObs[i].Sys != sys_t::GPS) continue; //只有GPS
        //if (obs->SatObs[i].Sys != sys_t::BDS) continue; //只有BDS
        double xsr = pos->Pos[0] - obs->SatPVT[i].Pos[0];
        double ysr = pos->Pos[1] - obs->SatPVT[i].Pos[1];
        double zsr = pos->Pos[2] - obs->SatPVT[i].Pos[2];
        double ρ = sqrt(xsr * xsr + ysr * ysr + zsr * zsr);
        double ρdot = -(xsr * obs->SatPVT[i].V[0] + ysr * obs->SatPVT[i].V[1] + zsr * obs->SatPVT[i].V[2]) / ρ;
        B(i, 0) = xsr / ρ;
        B(i, 1) = ysr / ρ;
        B(i, 2) = zsr / ρ;
        B(i, 3) = 1;
        w(i, 0) = obs->SatObs[i].d[0] - (ρdot - CLIGHT * obs->SatPVT[i].ClkSft);
        valid_satnum++;
    }
    // 可用卫星数不足以解算，则退出函数
    if (valid_satnum < 4) return;

    Matrix N = B.transpose() * B;
    Matrix W = B.transpose() * w;
    // 最小二乘求解 
    Matrix N_inv = N.inverse();
    Matrix X = N_inv * W;
    // 存储测速结果
    pos->Vel[0] = X(0, 0);
    pos->Vel[1] = X(1, 0);
    pos->Vel[2] = X(2, 0);
    // 计算验后单位权中误差
    Matrix v = B * X - w;
    Matrix vTv = v.transpose() * v;
    pos->SigmaVel = sqrt(vTv(0, 0) / (valid_satnum - 4));
}

/*
*******************************************
函数名：结果输出函数
参数：pos        解算结果
      refpos     参考位置
      denu       定位误差
函数功能：将结果进行形式转换并打印在控制台
*******************************************
*/
void OutputResult(const epochPOS_t* pos, XYZ refpos, Matrix& denu)
{
    std::cout << pos->Time.Week << " " << std::setw(7) << std::setfill(' ') << setprecision(0) << pos->Time.SecOfWeek << " 卫星数:" << pos->SatNum << "   解算位置:" << std::fixed << std::setprecision(3)
        << std::setw(14) << std::setfill(' ') << pos->Pos[0] << " " << std::setw(13) << std::setfill(' ') << pos->Pos[1] << " "
        << std::setw(13) << std::setfill(' ') << pos->Pos[2] << "  PDOP:" << std::setw(6) << std::setfill(' ') << pos->PDOP << "  sigmaPos:"
        << std::setw(6) << std::setfill(' ') << pos->SigmaPos << "   解算速度: " << std::setw(6) << std::setfill(' ') << pos->Vel[0] << " "
        << std::setw(6) << std::setfill(' ') << pos->Vel[1] << " " << std::setw(6) << std::setfill(' ') << pos->Vel[2] << "  sigmaV:"
        << std::setw(5) << std::setfill(' ') << pos->SigmaVel << endl;

   /* // 计算定位误差
    BLH blh;
    XYZToBLH(refpos, &blh, R_CGS2K, F_CGS2K);
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

    Matrix dxyz(3, 1, 0);
    dxyz(0, 0) = pos->Pos[0] - refpos.x;
    dxyz(1, 0) = pos->Pos[1] - refpos.y;
    dxyz(2, 0) = pos->Pos[2] - refpos.z;
    denu = Mat * dxyz;
    std::cout << "  定位误差(ENU):" << std::fixed << std::setprecision(3) << std::setw(5) << std::setfill(' ') << denu(0, 0) << " "
        << std::setw(5) << std::setfill(' ') << denu(1, 0) << " "
        << std::setw(5) << std::setfill(' ') << denu(2, 0) << std::endl << std::endl;
        */
}