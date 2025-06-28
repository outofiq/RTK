#include "header.h"

/*
*********************************************************************
函数名：计算GPS中间量（卫星的位置、速度，钟差、钟速，硬件延迟）的函数
参数：Prn     卫星号
	  t       需知中间量对应的时刻（GPS时）
	  eph     GPS广播星历
	  pos     所得中间量结果存储体
返回值：true-计算成功 false-计算失败
函数功能：根据卫星广播星历来计算中间量
*********************************************************************
*/
bool CalGPSPos(const int Prn,const GPSTIME* t,const gpseph_t* eph, SatPos_t* pos)
{
    double delta_t = diffTime(t, &eph->TOE);
	// 星历是否过期或健康判断
	if (eph->health == 1 || (fabs(delta_t) > 2*SECPERHOUR))
	{
		pos->Valid = false;
		return false;
	}
    else pos->Valid = true;

	// 计算卫星位置
	// 计算轨道长半轴
	double A = eph->SqrtA * eph->SqrtA;
	// 计算平均运动角速度
	double n0 = sqrt(GM_WGS / (A * A * A));// rad/s
	// 计算相对于星历参考历元的时间
	double tk = delta_t;
	// 对平均运动角速度进行改正
	double n = n0 + eph->DeltaN;// rad/s
	// 计算平近点角
	double Mk = eph->M0 + n * tk;// rad
	// 计算偏近点角（迭代）
	double Ek = 0;// rad 
	double Et = Mk;// 迭代器
	while (fabs(Et - Ek) > 1e-12)
	{
		Et = Ek;
		Ek = Mk + eph->ecc * sin(Ek);
	}
	// 计算真近点角
	double vk = atan2(sqrt(1 - eph->ecc * eph->ecc) * sin(Ek), cos(Ek) - eph->ecc);
	// 计算升交角距
	double PHIk = vk + eph->omega;
	// 计算二阶调和改正数
	double duk = eph->Cus * sin(2 * PHIk) + eph->Cuc * cos(2 * PHIk);// 升交角距改正数
	double drk = eph->Crs * sin(2 * PHIk) + eph->Crc * cos(2 * PHIk);// 向径改正数
	double dik = eph->Cis * sin(2 * PHIk) + eph->Cic * cos(2 * PHIk);// 轨道倾角改正数
	// 计算经过改正的升交角距
	double uk = PHIk + duk;
	// 计算经过改正的向径
	double rk = A * (1 - eph->ecc * cos(Ek)) + drk;
	// 计算经过改正的轨道倾角
	double ik = eph->i0 + dik + eph->iDot * tk;
	// 计算卫星在轨道平面上的位置
	double xk1 = rk * cos(uk);
	double yk1 = rk * sin(uk);
	// 计算改正后的升交点经度
	double OMEGAk = eph->OMEGA0 + (eph->OMEGADot - Omega_WGS) * tk - Omega_WGS * eph->TOE.SecOfWeek;
	// 计算在地固坐标系下的位置
	double xk = xk1 * cos(OMEGAk) - yk1 * cos(ik) * sin(OMEGAk);
	double yk = xk1 * sin(OMEGAk) + yk1 * cos(ik) * cos(OMEGAk);
	double zk = yk1 * sin(ik);
	pos->Pos[0] = xk;
	pos->Pos[1] = yk;
	pos->Pos[2] = zk;

	// 计算卫星运动速度
	double Ekdot = n / (1 - eph->ecc * cos(Ek));
	double PHIkdot = sqrt((1 + eph->ecc) / (1 - eph->ecc)) * (cos(vk / 2) * cos(vk / 2) / (cos(Ek / 2) * cos(Ek / 2))) * Ekdot;
	double ukdot = 2 * (eph->Cus * cos(2 * PHIk) - eph->Cuc * sin(2 * PHIk)) * PHIkdot + PHIkdot;
	double rkdot = A * eph->ecc * sin(Ek) * Ekdot + 2 * (eph->Crs * cos(2 * PHIk) - eph->Crc * sin(2 * PHIk)) * PHIkdot;
	double ikdot = eph->iDot + 2 * (eph->Cis * cos(2 * PHIk) - eph->Cic * sin(2 * PHIk)) * PHIkdot;
	double OMEGAkdot = eph->OMEGADot - Omega_WGS;
	// 列出Rdot矩阵
	Matrix Rdot(3, 4, 0);
	Rdot(0, 0) = cos(OMEGAk);
	Rdot(0, 1) = -sin(OMEGAk) * cos(ik);
	Rdot(0, 2) = -(xk1 * sin(OMEGAk) + yk1 * cos(OMEGAk) * cos(ik));
	Rdot(0, 3) = yk1 * sin(OMEGAk) * sin(ik);
	Rdot(1, 0) = sin(OMEGAk);
	Rdot(1, 1) = cos(OMEGAk) * cos(ik);
	Rdot(1, 2) = (xk1 * cos(OMEGAk) - yk1 * sin(OMEGAk) * cos(ik));
    Rdot(1, 3) = yk1 * cos(OMEGAk) * sin(ik);
	Rdot(2, 0) = 0.0;
	Rdot(2, 1) = sin(ik);
	Rdot(2, 2) = 0.0;
	Rdot(2, 3) = yk1 * cos(ik);
	// 列出辅助量
	double xk1dot = rkdot * cos(uk) - rk * ukdot * sin(uk);
	double yk1dot = rkdot * sin(uk) + rk * ukdot * cos(uk);
	Matrix xyomegai(4, 1, 0);
	xyomegai(0, 0) = xk1dot;
	xyomegai(1, 0) = yk1dot;
	xyomegai(2, 0) = OMEGAkdot;
	xyomegai(3, 0) = ikdot;
	Matrix V = Rdot * xyomegai;
	pos->V[0] = V(0, 0);
	pos->V[1] = V(1, 0);
	pos->V[2] = V(2, 0);

	// 计算钟差
	double dtr = -4.442807633e-10 * eph->ecc * eph->SqrtA * sin(Ek); // 相对论效应改正
    double t_toc = diffTime(t, &eph->TOC);
    double dtsv = eph->ClkBias + eph->ClkDrift * t_toc + eph->ClkDriftRate * t_toc * t_toc + dtr;
	pos->ClkOft = dtsv;

	// 计算钟速
	double dtrdot = -4.442807633e-10 * eph->ecc * eph->SqrtA * cos(Ek) * Ekdot; // 相对论效应改正
	double dtsvdot = eph->ClkDrift + 2 * eph->ClkDriftRate * t_toc + dtrdot;
	pos->ClkSft = dtsvdot;

	// 硬件延迟项
	pos->Tgd1 = eph->TGD1;
	pos->Tgd2 = eph->TGD2;

	return true;
}

/*
*********************************************************************
函数名：计算BDS中间量（卫星的位置、速度，钟差、钟速，硬件延迟）的函数
参数：Prn     卫星号
	  t       需知中间量对应的时刻（GPS时）
	  eph     BDS广播星历
	  pos     所得中间量结果存储体
返回值：true-计算成功 false-计算失败
函数功能：根据卫星广播星历来计算中间量
*********************************************************************
*/
bool CalBDSPos(const int Prn,const GPSTIME* t,const gpseph_t* eph, SatPos_t* pos)
{
    // 把t转化到BDS时下
    GPSTIME tb;
    tb.Week = t->Week - 1356;
    tb.SecOfWeek = t->SecOfWeek - 14;
    double delta_t = diffTime(&tb, &eph->TOE);
    // 星历是否过期或健康判断
    if (eph->health == 1 || (fabs(delta_t) > SECPERHOUR))
    {
        pos->Valid = false;
        return false;
    }
    else pos->Valid = true;

    // 计算卫星位置
    // 计算轨道长半轴
    double A = eph->SqrtA * eph->SqrtA;
    // 计算平均运动角速度
    double n0 = sqrt(GM_BDS / (A * A * A));// rad/s
    // 计算相对于星历参考历元的时间
    double tk = delta_t;
    // 对平均运动角速度进行改正
    double n = n0 + eph->DeltaN;// rad/s
    // 计算平近点角
    double Mk = eph->M0 + n * tk;// rad
    // 计算偏近点角（迭代）
    double Ek = 0;// rad 
    double Et = Mk;//迭代器
    while (fabs(Et - Ek) > 1e-12)
    {
        Et = Ek;
        Ek = Mk + eph->ecc * sin(Ek);
    }
    // 计算真近点角
    double vk = atan2(sqrt(1 - eph->ecc * eph->ecc) * sin(Ek), (cos(Ek) - eph->ecc));
    // 计算升交角距
    double PHIk = vk + eph->omega;
    // 计算二阶调和改正数
    double duk = eph->Cus * sin(2 * PHIk) + eph->Cuc * cos(2 * PHIk);//升交角距改正数
    double drk = eph->Crs * sin(2 * PHIk) + eph->Crc * cos(2 * PHIk);//向径改正数
    double dik = eph->Cis * sin(2 * PHIk) + eph->Cic * cos(2 * PHIk);//轨道倾角改正数
    // 计算经过改正的升交角距
    double uk = PHIk + duk;
    // 计算经过改正的向径
    double rk = A * (1 - eph->ecc * cos(Ek)) + drk;
    // 计算经过改正的轨道倾角
    double ik = eph->i0 + dik + eph->iDot * tk;
    // 计算卫星在轨道平面上的位置
    double xk1 = rk * cos(uk);
    double yk1 = rk * sin(uk);
    double OMEGAk;

    // 根据卫星种类不同进行选择
    // GEO卫星
    if ((Prn >= 1 && Prn <= 5) || (Prn >= 59 && Prn <= 63))
    {
        // 计算历元升交点经度
        OMEGAk = eph->OMEGA0 + eph->OMEGADot * tk - Omega_BDS * eph->TOE.SecOfWeek;
        // 卫星在自定义的坐标系中的坐标
        double xgk = xk1 * cos(OMEGAk) - yk1 * cos(ik) * sin(OMEGAk);
        double ygk = xk1 * sin(OMEGAk) + yk1 * cos(ik) * cos(OMEGAk);
        double zgk = yk1 * sin(ik);
        Matrix gk(3, 1, 0);
        gk(0, 0) = xgk;
        gk(1, 0) = ygk;
        gk(2, 0) = zgk;
        // 投影矩阵
        Matrix Rx(3, 3, 0);
        Rx(0, 0) = 1;
        Rx(0, 1) = 0;
        Rx(0, 2) = 0;
        Rx(1, 0) = 0;
        Rx(1, 1) = cos(-5 * D2R);
        Rx(1, 2) = sin(-5 * D2R);
        Rx(2, 0) = 0;
        Rx(2, 1) = -sin(-5 * D2R);
        Rx(2, 2) = cos(-5 * D2R);

        Matrix Rz(3, 3, 0);
        Rz(0, 0) = cos(Omega_BDS * tk);
        Rz(0, 1) = sin(Omega_BDS * tk);
        Rz(0, 2) = 0;
        Rz(1, 0) = -sin(Omega_BDS * tk);
        Rz(1, 1) = cos(Omega_BDS * tk);
        Rz(1, 2) = 0;
        Rz(2, 0) = 0;
        Rz(2, 1) = 0;
        Rz(2, 2) = 1;
       // GEO卫星在BDCS坐标系中的坐标
        Matrix posK = Rz * Rx * gk;
        pos->Pos[0] = posK(0, 0);
        pos->Pos[1] = posK(1, 0);
        pos->Pos[2] = posK(2, 0);

        // 计算卫星运动速度
        double Ekdot = n / (1 - eph->ecc * cos(Ek));
        double PHIkdot = sqrt((1 + eph->ecc) / (1 - eph->ecc)) * (cos(vk / 2) * cos(vk / 2) / (cos(Ek / 2) * cos(Ek / 2))) * Ekdot;
        double ukdot = 2 * (eph->Cus * cos(2 * PHIk) - eph->Cuc * sin(2 * PHIk)) * PHIkdot + PHIkdot;
        double rkdot = A * eph->ecc * sin(Ek) * Ekdot + 2 * (eph->Crs * cos(2 * PHIk) - eph->Crc * sin(2 * PHIk)) * PHIkdot;
        double ikdot = eph->iDot + 2 * (eph->Cis * cos(2 * PHIk) - eph->Cic * sin(2 * PHIk)) * PHIkdot;
        // 计算卫星在轨道平面上的速度
        double xk1dot = rkdot * cos(uk) - rk * ukdot * sin(uk);
        double yk1dot = rkdot * sin(uk) + rk * ukdot * cos(uk);

        double OMEGAkdot = eph->OMEGADot;
        //计算卫星在自定义旋转坐标系中的速度
        double xgkdot = -ygk * OMEGAkdot - (yk1dot * cos(ik) - zgk * ikdot) * sin(OMEGAk) + xk1dot * cos(OMEGAk);
        double ygkdot = xgk * OMEGAkdot + (yk1dot * cos(ik) - zgk * ikdot) * cos(OMEGAk) + xk1dot * sin(OMEGAk);
        double zgkdot = yk1dot * sin(ik) + yk1 * ikdot * cos(ik);

        Matrix gkdot(3, 1, 0);
        gkdot(0, 0) = xgkdot;
        gkdot(1, 0) = ygkdot;
        gkdot(2, 0) = zgkdot;

        Matrix Rzdot(3, 3, 0);
        Rzdot(0, 0) = -sin(Omega_BDS * tk) * Omega_BDS;
        Rzdot(0, 1) = cos(Omega_BDS * tk) * Omega_BDS;
        Rzdot(0, 2) = 0;
        Rzdot(1, 0) = -cos(Omega_BDS * tk) * Omega_BDS;
        Rzdot(1, 1) = -sin(Omega_BDS * tk) * Omega_BDS;
        Rzdot(1, 2) = 0;
        Rzdot(2, 0) = 0;
        Rzdot(2, 1) = 0;
        Rzdot(2, 2) = 0;
        // 计算GEO卫星在BDCS坐标系下的速度
        Matrix Vk = Rzdot * Rx * gk + Rz * Rx * gkdot;
        pos->V[0] = Vk(0, 0);
        pos->V[1] = Vk(1, 0);
        pos->V[2] = Vk(2, 0);

        // 计算钟差
        double dtr = -4.442807633e-10 * eph->ecc * eph->SqrtA * sin(Ek);// 相对论效应改正
        double t_toc = diffTime(&tb, &eph->TOC);
        double dtsv = eph->ClkBias + eph->ClkDrift * t_toc + eph->ClkDriftRate * t_toc * t_toc + dtr;
        pos->ClkOft = dtsv;

        // 计算钟速
        double dtrdot = -4.442807633e-10 * eph->ecc * eph->SqrtA * cos(Ek) * Ekdot;// 相对论效应改正
        double ddsv = eph->ClkDrift + 2 * eph->ClkDriftRate * t_toc + dtrdot;
        pos->ClkSft = ddsv;
    }
    // MEO或IGSO卫星
    else if (Prn > 5 && Prn < 59)
    {
        // 计算改正后的升交点经度
        OMEGAk = eph->OMEGA0 + (eph->OMEGADot - Omega_BDS) * tk - Omega_BDS * eph->TOE.SecOfWeek;
        // 计算在地固坐标系下的位置
        double xk = xk1 * cos(OMEGAk) - yk1 * cos(ik) * sin(OMEGAk);
        double yk = xk1 * sin(OMEGAk) + yk1 * cos(ik) * cos(OMEGAk);
        double zk = yk1 * sin(ik);
        pos->Pos[0] = xk;
        pos->Pos[1] = yk;
        pos->Pos[2] = zk;

        // 计算卫星运动速度
        double Ekdot = n / (1 - eph->ecc * cos(Ek));
        double PHIkdot = sqrt((1 + eph->ecc) / (1 - eph->ecc)) * (cos(vk / 2) * cos(vk / 2) / (cos(Ek / 2) * cos(Ek / 2))) * Ekdot;
        double ukdot = 2 * (eph->Cus * cos(2 * PHIk) - eph->Cuc * sin(2 * PHIk)) * PHIkdot + PHIkdot;
        double rkdot = A * eph->ecc * sin(Ek) * Ekdot + 2 * (eph->Crs * cos(2 * PHIk) - eph->Crc * sin(2 * PHIk)) * PHIkdot;
        double ikdot = eph->iDot + 2 * (eph->Cis * cos(2 * PHIk) - eph->Cic * sin(2 * PHIk)) * PHIkdot;
        // 计算卫星在轨道平面上的速度
        double xk1dot = rkdot * cos(uk) - rk * ukdot * sin(uk);
        double yk1dot = rkdot * sin(uk) + rk * ukdot * cos(uk);
        //计算升交点经度变化率（地固系）
        double OMEGAkdot = eph->OMEGADot - Omega_BDS;
        //计算MEO/IGSO卫星在BDCS坐标系中的速度
        double xkdot = -yk * OMEGAkdot - (yk1dot * cos(ik) - zk * ikdot) * sin(OMEGAk) + xk1dot * cos(OMEGAk);
        double ykdot = xk * OMEGAkdot + (yk1dot * cos(ik) - zk * ikdot) * cos(OMEGAk) + xk1dot * sin(OMEGAk);
        double zkdot = yk1dot * sin(ik) + yk1 * ikdot * cos(ik);
        pos->V[0] = xkdot;
        pos->V[1] = ykdot;
        pos->V[2] = zkdot;

        // 计算钟差
        double dtr = -4.442807633e-10 * eph->ecc * eph->SqrtA * sin(Ek);// 相对论效应改正
        double t_toc = diffTime(&tb, &eph->TOC);
        double dtsv = eph->ClkBias + eph->ClkDrift * t_toc + eph->ClkDriftRate * t_toc * t_toc + dtr;
        pos->ClkOft = dtsv;

        // 计算钟速
        double dtrdot = -4.442807633e-10 * eph->ecc * eph->SqrtA * cos(Ek) * Ekdot;// 相对论效应改正
        double ddsv = eph->ClkDrift + 2 * eph->ClkDriftRate * t_toc + dtrdot;
        pos->ClkSft = ddsv;
    }
    // 硬件延迟项
    pos->Tgd1 = eph->TGD1;
    pos->Tgd2 = eph->TGD2;

    return true;
}

/*
*********************************************************************
函数名：计算信号发射时刻的卫星位置、速度、钟差、钟速
参数：obs     当前历元观测数据
      GPSEph  GPS星历
      BDSEph  BDS星历
      xyz     接收机的XYZ坐标(m)
函数功能：计算信号发射时刻的卫星位置、速度、钟差、钟速，并改正位置和速度，计算高度角、方位角、对流层延迟
*********************************************************************
*/
void CalSatPos(epoch_t* obs, gpseph_t* GPSEph, gpseph_t* BDSEph, XYZ xyz)
{
    GPSTIME t_clock = obs->Time;
    for (int i = 0; i < obs->SatNum; i++)
    {
        GPSTIME t_AtSignalTrans;
        double dt = 0;// 迭代计算卫星钟差阈值判断
        if (obs->SatObs[i].Sys == sys_t::GPS)
        {
            gpseph_t* eph = GPSEph + obs->SatObs[i].Prn - 1;// 取用本组观测值对应的卫星的星历
            obs->SatPVT[i].ClkOft = 0;// 初始钟差为0
            do
            {   // 计算卫星信号发射时刻
                t_AtSignalTrans.Week = t_clock.Week;
                t_AtSignalTrans.SecOfWeek = t_clock.SecOfWeek - obs->SatObs[i].p[0] / CLIGHT- obs->SatPVT[i].ClkOft;
                // 计算卫星钟差
                double t_toc = diffTime(&t_AtSignalTrans, &eph->TOC);
                double st_tmp = eph->ClkBias + eph->ClkDrift * t_toc + eph->ClkDriftRate * t_toc * t_toc;
                // 更新判断值
                dt = fabs(st_tmp - obs->SatPVT[i].ClkOft);
                // 更新卫星钟差
                obs->SatPVT[i].ClkOft = st_tmp;
            } while (dt > 1e-12);
            // 卫星钟差迭代计算完毕
            // 计算卫星位置、速度、钟差和钟速
            if (CalGPSPos(obs->SatObs[i].Prn, &t_AtSignalTrans, eph, obs->SatPVT + i))
            {
                // 计算信号传播时间
                double t_trans = sqrt((obs->SatPVT[i].Pos[0] - xyz.x) * (obs->SatPVT[i].Pos[0] - xyz.x) + (obs->SatPVT[i].Pos[1] - xyz.y) * (obs->SatPVT[i].Pos[1] - xyz.y) + (obs->SatPVT[i].Pos[2] - xyz.z) * (obs->SatPVT[i].Pos[2] - xyz.z)) / CLIGHT;
                // 地球自转改正
                double alpha = Omega_WGS * t_trans;

                Matrix Rz(3, 3, 0);
                Rz(0, 0) = cos(alpha);
                Rz(0, 1) = sin(alpha);
                Rz(0, 2) = 0;
                Rz(1, 0) = -sin(alpha);
                Rz(1, 1) = cos(alpha);
                Rz(1, 2) = 0;
                Rz(2, 0) = 0;
                Rz(2, 1) = 0;
                Rz(2, 2) = 1;

                Matrix pos(3, 1, 0);
                pos(0, 0) = obs->SatPVT[i].Pos[0];
                pos(1, 0) = obs->SatPVT[i].Pos[1];
                pos(2, 0) = obs->SatPVT[i].Pos[2];

                Matrix vel(3, 1, 0);
                vel(0, 0) = obs->SatPVT[i].V[0];
                vel(1, 0) = obs->SatPVT[i].V[1];
                vel(2, 0) = obs->SatPVT[i].V[2];

                // 位置改正
                Matrix pos_new = Rz * pos;
                obs->SatPVT[i].Pos[0] = pos_new(0, 0);
                obs->SatPVT[i].Pos[1] = pos_new(1, 0);
                obs->SatPVT[i].Pos[2] = pos_new(2, 0);
                // 速度改正
                Matrix vel_new = Rz * vel;
                obs->SatPVT[i].V[0] = vel_new(0, 0);
                obs->SatPVT[i].V[1] = vel_new(1, 0);
                obs->SatPVT[i].V[2] = vel_new(2, 0);

                // 计算卫星的高度角和方位角
                BLH blh;
                XYZToBLH(xyz, &blh, R_WGS84, F_WGS84);
                // 计算测站地平坐标转换矩阵
                double B = blh.B;
                double L = blh.L;
                double H = blh.H;
                double sinL = sin(L);
                double cosL = cos(L);
                double sinB = sin(B);
                double cosB = cos(B);
                Matrix Mat(3, 3, 0);
                Mat(0, 0) = -sinL;
                Mat(0, 1) = cosL;
                Mat(0, 2) = 0;
                Mat(1, 0) = -sinB * cosL;
                Mat(1, 1) = -sinB * sinL;
                Mat(1, 2) = cosB;
                Mat(2, 0) = cosB * cosL;
                Mat(2, 1) = cosB * sinL;
                Mat(2, 2) = sinB;

                Matrix dxyz(3, 1, 0);
                dxyz(0, 0) = obs->SatPVT[i].Pos[0] - xyz.x;
                dxyz(1, 0) = obs->SatPVT[i].Pos[1] - xyz.y;
                dxyz(2, 0) = obs->SatPVT[i].Pos[2] - xyz.z;

                Matrix dENU = Mat * dxyz;
                double len = sqrt(dENU(0, 0) * dENU(0, 0) + dENU(1, 0) * dENU(1, 0) + dENU(2, 0) * dENU(2, 0));
                obs->SatPVT[i].Elevation = asin(dENU(2, 0) / len);
                obs->SatPVT[i].Azimuth = atan2(dENU(0, 0), dENU(1, 0));

                // 计算对流层改正
                obs->SatPVT[i].TropCorr = hopfield(blh.H, obs->SatPVT[i].Elevation);
            }
        }
        else if (obs->SatObs[i].Sys == sys_t::BDS)
        {

            gpseph_t* eph = BDSEph + obs->SatObs[i].Prn - 1;// 取用本组观测值对应的卫星的星历
            obs->SatPVT[i].ClkOft = 0;// 初始钟差为0
            do
            {   // 计算卫星信号发射时刻
                t_AtSignalTrans.Week = t_clock.Week;
                t_AtSignalTrans.SecOfWeek = t_clock.SecOfWeek - obs->SatObs[i].p[0] / CLIGHT - obs->SatPVT[i].ClkOft;
                // 计算卫星钟差
                GPSTIME tb;
                tb.Week = t_AtSignalTrans.Week - 1356;
                tb.SecOfWeek = t_AtSignalTrans.SecOfWeek - 14;
                double t_toc = diffTime(&tb, &eph->TOC);
                double st_tmp = eph->ClkBias + eph->ClkDrift * t_toc + eph->ClkDriftRate * t_toc * t_toc;
                // 更新判断值
                dt = fabs(st_tmp - obs->SatPVT[i].ClkOft);
                // 更新卫星钟差
                obs->SatPVT[i].ClkOft = st_tmp;
            } while (dt > 1e-12);
            // 卫星钟差迭代计算完毕
            // 计算卫星位置、速度、钟差和钟速
            if (CalBDSPos(obs->SatObs[i].Prn, &t_AtSignalTrans, eph, obs->SatPVT + i))
            {
                // 计算信号传播时间
                double t_trans = sqrt((obs->SatPVT[i].Pos[0] - xyz.x) * (obs->SatPVT[i].Pos[0] - xyz.x) + (obs->SatPVT[i].Pos[1] - xyz.y) * (obs->SatPVT[i].Pos[1] - xyz.y) + (obs->SatPVT[i].Pos[2] - xyz.z) * (obs->SatPVT[i].Pos[2] - xyz.z)) / CLIGHT;
                // 地球自转改正
                double alpha = Omega_BDS * t_trans;
                Matrix Rz(3, 3, 0);
                Rz(0, 0) = cos(alpha);
                Rz(0, 1) = sin(alpha);
                Rz(0, 2) = 0;
                Rz(1, 0) = -sin(alpha);
                Rz(1, 1) = cos(alpha);
                Rz(1, 2) = 0;
                Rz(2, 0) = 0;
                Rz(2, 1) = 0;
                Rz(2, 2) = 1;

                Matrix pos(3, 1, 0);
                pos(0, 0) = obs->SatPVT[i].Pos[0];
                pos(1, 0) = obs->SatPVT[i].Pos[1];
                pos(2, 0) = obs->SatPVT[i].Pos[2];

                Matrix vel(3, 1, 0);
                vel(0, 0) = obs->SatPVT[i].V[0];
                vel(1, 0) = obs->SatPVT[i].V[1];
                vel(2, 0) = obs->SatPVT[i].V[2];

                // 位置改正
                Matrix pos_new = Rz * pos;
                obs->SatPVT[i].Pos[0] = pos_new(0, 0);
                obs->SatPVT[i].Pos[1] = pos_new(1, 0);
                obs->SatPVT[i].Pos[2] = pos_new(2, 0);
                // 速度改正
                Matrix vel_new = Rz * vel;
                obs->SatPVT[i].V[0] = vel_new(0, 0);
                obs->SatPVT[i].V[1] = vel_new(1, 0);
                obs->SatPVT[i].V[2] = vel_new(2, 0);

                // 计算卫星的高度角和方位角
                BLH blh;
                XYZToBLH(xyz, &blh, R_CGS2K, F_CGS2K);
                // 计算测站地平坐标转换矩阵
                double B = blh.B;
                double L = blh.L;
                double H = blh.H;
                Matrix Mat(3, 3, 0);
                Mat(0, 0) = -sin(L);
                Mat(0, 1) = cos(L);
                Mat(0, 2) = 0;
                Mat(1, 0) = -sin(B) * cos(L);
                Mat(1, 1) = -sin(B) * sin(L);
                Mat(1, 2) = cos(B);
                Mat(2, 0) = cos(B) * cos(L);
                Mat(2, 1) = cos(B) * sin(L);
                Mat(2, 2) = sin(B);

                Matrix dxyz(3, 1, 0);
                dxyz(0, 0) = obs->SatPVT[i].Pos[0] - xyz.x;
                dxyz(1, 0) = obs->SatPVT[i].Pos[1] - xyz.y;
                dxyz(2, 0) = obs->SatPVT[i].Pos[2] - xyz.z;

                Matrix dENU = Mat * dxyz;
                double len = sqrt(dENU(0, 0) * dENU(0, 0) + dENU(1, 0) * dENU(1, 0) + dENU(2, 0) * dENU(2, 0));
                obs->SatPVT[i].Elevation = asin(dENU(2, 0) / len);
                obs->SatPVT[i].Azimuth = atan2(dENU(0, 0), dENU(1, 0));

                // 计算对流层改正
                obs->SatPVT[i].TropCorr = hopfield(blh.H, obs->SatPVT[i].Elevation);
            }
        }
    }
}