#include "header.h"

/*
*********************************************************************
������������GPS�м��������ǵ�λ�á��ٶȣ��Ӳ���٣�Ӳ���ӳ٣��ĺ���
������Prn     ���Ǻ�
	  t       ��֪�м�����Ӧ��ʱ�̣�GPSʱ��
	  eph     GPS�㲥����
	  pos     �����м�������洢��
����ֵ��true-����ɹ� false-����ʧ��
�������ܣ��������ǹ㲥�����������м���
*********************************************************************
*/
bool CalGPSPos(const int Prn,const GPSTIME* t,const gpseph_t* eph, SatPos_t* pos)
{
    double delta_t = diffTime(t, &eph->TOE);
	// �����Ƿ���ڻ򽡿��ж�
	if (eph->health == 1 || (fabs(delta_t) > 2*SECPERHOUR))
	{
		pos->Valid = false;
		return false;
	}
    else pos->Valid = true;

	// ��������λ��
	// ������������
	double A = eph->SqrtA * eph->SqrtA;
	// ����ƽ���˶����ٶ�
	double n0 = sqrt(GM_WGS / (A * A * A));// rad/s
	// ��������������ο���Ԫ��ʱ��
	double tk = delta_t;
	// ��ƽ���˶����ٶȽ��и���
	double n = n0 + eph->DeltaN;// rad/s
	// ����ƽ�����
	double Mk = eph->M0 + n * tk;// rad
	// ����ƫ����ǣ�������
	double Ek = 0;// rad 
	double Et = Mk;// ������
	while (fabs(Et - Ek) > 1e-12)
	{
		Et = Ek;
		Ek = Mk + eph->ecc * sin(Ek);
	}
	// ����������
	double vk = atan2(sqrt(1 - eph->ecc * eph->ecc) * sin(Ek), cos(Ek) - eph->ecc);
	// ���������Ǿ�
	double PHIk = vk + eph->omega;
	// ������׵��͸�����
	double duk = eph->Cus * sin(2 * PHIk) + eph->Cuc * cos(2 * PHIk);// �����Ǿ������
	double drk = eph->Crs * sin(2 * PHIk) + eph->Crc * cos(2 * PHIk);// �򾶸�����
	double dik = eph->Cis * sin(2 * PHIk) + eph->Cic * cos(2 * PHIk);// �����Ǹ�����
	// ���㾭�������������Ǿ�
	double uk = PHIk + duk;
	// ���㾭����������
	double rk = A * (1 - eph->ecc * cos(Ek)) + drk;
	// ���㾭�������Ĺ�����
	double ik = eph->i0 + dik + eph->iDot * tk;
	// ���������ڹ��ƽ���ϵ�λ��
	double xk1 = rk * cos(uk);
	double yk1 = rk * sin(uk);
	// ���������������㾭��
	double OMEGAk = eph->OMEGA0 + (eph->OMEGADot - Omega_WGS) * tk - Omega_WGS * eph->TOE.SecOfWeek;
	// �����ڵع�����ϵ�µ�λ��
	double xk = xk1 * cos(OMEGAk) - yk1 * cos(ik) * sin(OMEGAk);
	double yk = xk1 * sin(OMEGAk) + yk1 * cos(ik) * cos(OMEGAk);
	double zk = yk1 * sin(ik);
	pos->Pos[0] = xk;
	pos->Pos[1] = yk;
	pos->Pos[2] = zk;

	// ���������˶��ٶ�
	double Ekdot = n / (1 - eph->ecc * cos(Ek));
	double PHIkdot = sqrt((1 + eph->ecc) / (1 - eph->ecc)) * (cos(vk / 2) * cos(vk / 2) / (cos(Ek / 2) * cos(Ek / 2))) * Ekdot;
	double ukdot = 2 * (eph->Cus * cos(2 * PHIk) - eph->Cuc * sin(2 * PHIk)) * PHIkdot + PHIkdot;
	double rkdot = A * eph->ecc * sin(Ek) * Ekdot + 2 * (eph->Crs * cos(2 * PHIk) - eph->Crc * sin(2 * PHIk)) * PHIkdot;
	double ikdot = eph->iDot + 2 * (eph->Cis * cos(2 * PHIk) - eph->Cic * sin(2 * PHIk)) * PHIkdot;
	double OMEGAkdot = eph->OMEGADot - Omega_WGS;
	// �г�Rdot����
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
	// �г�������
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

	// �����Ӳ�
	double dtr = -4.442807633e-10 * eph->ecc * eph->SqrtA * sin(Ek); // �����ЧӦ����
    double t_toc = diffTime(t, &eph->TOC);
    double dtsv = eph->ClkBias + eph->ClkDrift * t_toc + eph->ClkDriftRate * t_toc * t_toc + dtr;
	pos->ClkOft = dtsv;

	// ��������
	double dtrdot = -4.442807633e-10 * eph->ecc * eph->SqrtA * cos(Ek) * Ekdot; // �����ЧӦ����
	double dtsvdot = eph->ClkDrift + 2 * eph->ClkDriftRate * t_toc + dtrdot;
	pos->ClkSft = dtsvdot;

	// Ӳ���ӳ���
	pos->Tgd1 = eph->TGD1;
	pos->Tgd2 = eph->TGD2;

	return true;
}

/*
*********************************************************************
������������BDS�м��������ǵ�λ�á��ٶȣ��Ӳ���٣�Ӳ���ӳ٣��ĺ���
������Prn     ���Ǻ�
	  t       ��֪�м�����Ӧ��ʱ�̣�GPSʱ��
	  eph     BDS�㲥����
	  pos     �����м�������洢��
����ֵ��true-����ɹ� false-����ʧ��
�������ܣ��������ǹ㲥�����������м���
*********************************************************************
*/
bool CalBDSPos(const int Prn,const GPSTIME* t,const gpseph_t* eph, SatPos_t* pos)
{
    // ��tת����BDSʱ��
    GPSTIME tb;
    tb.Week = t->Week - 1356;
    tb.SecOfWeek = t->SecOfWeek - 14;
    double delta_t = diffTime(&tb, &eph->TOE);
    // �����Ƿ���ڻ򽡿��ж�
    if (eph->health == 1 || (fabs(delta_t) > SECPERHOUR))
    {
        pos->Valid = false;
        return false;
    }
    else pos->Valid = true;

    // ��������λ��
    // ������������
    double A = eph->SqrtA * eph->SqrtA;
    // ����ƽ���˶����ٶ�
    double n0 = sqrt(GM_BDS / (A * A * A));// rad/s
    // ��������������ο���Ԫ��ʱ��
    double tk = delta_t;
    // ��ƽ���˶����ٶȽ��и���
    double n = n0 + eph->DeltaN;// rad/s
    // ����ƽ�����
    double Mk = eph->M0 + n * tk;// rad
    // ����ƫ����ǣ�������
    double Ek = 0;// rad 
    double Et = Mk;//������
    while (fabs(Et - Ek) > 1e-12)
    {
        Et = Ek;
        Ek = Mk + eph->ecc * sin(Ek);
    }
    // ����������
    double vk = atan2(sqrt(1 - eph->ecc * eph->ecc) * sin(Ek), (cos(Ek) - eph->ecc));
    // ���������Ǿ�
    double PHIk = vk + eph->omega;
    // ������׵��͸�����
    double duk = eph->Cus * sin(2 * PHIk) + eph->Cuc * cos(2 * PHIk);//�����Ǿ������
    double drk = eph->Crs * sin(2 * PHIk) + eph->Crc * cos(2 * PHIk);//�򾶸�����
    double dik = eph->Cis * sin(2 * PHIk) + eph->Cic * cos(2 * PHIk);//�����Ǹ�����
    // ���㾭�������������Ǿ�
    double uk = PHIk + duk;
    // ���㾭����������
    double rk = A * (1 - eph->ecc * cos(Ek)) + drk;
    // ���㾭�������Ĺ�����
    double ik = eph->i0 + dik + eph->iDot * tk;
    // ���������ڹ��ƽ���ϵ�λ��
    double xk1 = rk * cos(uk);
    double yk1 = rk * sin(uk);
    double OMEGAk;

    // �����������಻ͬ����ѡ��
    // GEO����
    if ((Prn >= 1 && Prn <= 5) || (Prn >= 59 && Prn <= 63))
    {
        // ������Ԫ�����㾭��
        OMEGAk = eph->OMEGA0 + eph->OMEGADot * tk - Omega_BDS * eph->TOE.SecOfWeek;
        // �������Զ��������ϵ�е�����
        double xgk = xk1 * cos(OMEGAk) - yk1 * cos(ik) * sin(OMEGAk);
        double ygk = xk1 * sin(OMEGAk) + yk1 * cos(ik) * cos(OMEGAk);
        double zgk = yk1 * sin(ik);
        Matrix gk(3, 1, 0);
        gk(0, 0) = xgk;
        gk(1, 0) = ygk;
        gk(2, 0) = zgk;
        // ͶӰ����
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
       // GEO������BDCS����ϵ�е�����
        Matrix posK = Rz * Rx * gk;
        pos->Pos[0] = posK(0, 0);
        pos->Pos[1] = posK(1, 0);
        pos->Pos[2] = posK(2, 0);

        // ���������˶��ٶ�
        double Ekdot = n / (1 - eph->ecc * cos(Ek));
        double PHIkdot = sqrt((1 + eph->ecc) / (1 - eph->ecc)) * (cos(vk / 2) * cos(vk / 2) / (cos(Ek / 2) * cos(Ek / 2))) * Ekdot;
        double ukdot = 2 * (eph->Cus * cos(2 * PHIk) - eph->Cuc * sin(2 * PHIk)) * PHIkdot + PHIkdot;
        double rkdot = A * eph->ecc * sin(Ek) * Ekdot + 2 * (eph->Crs * cos(2 * PHIk) - eph->Crc * sin(2 * PHIk)) * PHIkdot;
        double ikdot = eph->iDot + 2 * (eph->Cis * cos(2 * PHIk) - eph->Cic * sin(2 * PHIk)) * PHIkdot;
        // ���������ڹ��ƽ���ϵ��ٶ�
        double xk1dot = rkdot * cos(uk) - rk * ukdot * sin(uk);
        double yk1dot = rkdot * sin(uk) + rk * ukdot * cos(uk);

        double OMEGAkdot = eph->OMEGADot;
        //�����������Զ�����ת����ϵ�е��ٶ�
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
        // ����GEO������BDCS����ϵ�µ��ٶ�
        Matrix Vk = Rzdot * Rx * gk + Rz * Rx * gkdot;
        pos->V[0] = Vk(0, 0);
        pos->V[1] = Vk(1, 0);
        pos->V[2] = Vk(2, 0);

        // �����Ӳ�
        double dtr = -4.442807633e-10 * eph->ecc * eph->SqrtA * sin(Ek);// �����ЧӦ����
        double t_toc = diffTime(&tb, &eph->TOC);
        double dtsv = eph->ClkBias + eph->ClkDrift * t_toc + eph->ClkDriftRate * t_toc * t_toc + dtr;
        pos->ClkOft = dtsv;

        // ��������
        double dtrdot = -4.442807633e-10 * eph->ecc * eph->SqrtA * cos(Ek) * Ekdot;// �����ЧӦ����
        double ddsv = eph->ClkDrift + 2 * eph->ClkDriftRate * t_toc + dtrdot;
        pos->ClkSft = ddsv;
    }
    // MEO��IGSO����
    else if (Prn > 5 && Prn < 59)
    {
        // ���������������㾭��
        OMEGAk = eph->OMEGA0 + (eph->OMEGADot - Omega_BDS) * tk - Omega_BDS * eph->TOE.SecOfWeek;
        // �����ڵع�����ϵ�µ�λ��
        double xk = xk1 * cos(OMEGAk) - yk1 * cos(ik) * sin(OMEGAk);
        double yk = xk1 * sin(OMEGAk) + yk1 * cos(ik) * cos(OMEGAk);
        double zk = yk1 * sin(ik);
        pos->Pos[0] = xk;
        pos->Pos[1] = yk;
        pos->Pos[2] = zk;

        // ���������˶��ٶ�
        double Ekdot = n / (1 - eph->ecc * cos(Ek));
        double PHIkdot = sqrt((1 + eph->ecc) / (1 - eph->ecc)) * (cos(vk / 2) * cos(vk / 2) / (cos(Ek / 2) * cos(Ek / 2))) * Ekdot;
        double ukdot = 2 * (eph->Cus * cos(2 * PHIk) - eph->Cuc * sin(2 * PHIk)) * PHIkdot + PHIkdot;
        double rkdot = A * eph->ecc * sin(Ek) * Ekdot + 2 * (eph->Crs * cos(2 * PHIk) - eph->Crc * sin(2 * PHIk)) * PHIkdot;
        double ikdot = eph->iDot + 2 * (eph->Cis * cos(2 * PHIk) - eph->Cic * sin(2 * PHIk)) * PHIkdot;
        // ���������ڹ��ƽ���ϵ��ٶ�
        double xk1dot = rkdot * cos(uk) - rk * ukdot * sin(uk);
        double yk1dot = rkdot * sin(uk) + rk * ukdot * cos(uk);
        //���������㾭�ȱ仯�ʣ��ع�ϵ��
        double OMEGAkdot = eph->OMEGADot - Omega_BDS;
        //����MEO/IGSO������BDCS����ϵ�е��ٶ�
        double xkdot = -yk * OMEGAkdot - (yk1dot * cos(ik) - zk * ikdot) * sin(OMEGAk) + xk1dot * cos(OMEGAk);
        double ykdot = xk * OMEGAkdot + (yk1dot * cos(ik) - zk * ikdot) * cos(OMEGAk) + xk1dot * sin(OMEGAk);
        double zkdot = yk1dot * sin(ik) + yk1 * ikdot * cos(ik);
        pos->V[0] = xkdot;
        pos->V[1] = ykdot;
        pos->V[2] = zkdot;

        // �����Ӳ�
        double dtr = -4.442807633e-10 * eph->ecc * eph->SqrtA * sin(Ek);// �����ЧӦ����
        double t_toc = diffTime(&tb, &eph->TOC);
        double dtsv = eph->ClkBias + eph->ClkDrift * t_toc + eph->ClkDriftRate * t_toc * t_toc + dtr;
        pos->ClkOft = dtsv;

        // ��������
        double dtrdot = -4.442807633e-10 * eph->ecc * eph->SqrtA * cos(Ek) * Ekdot;// �����ЧӦ����
        double ddsv = eph->ClkDrift + 2 * eph->ClkDriftRate * t_toc + dtrdot;
        pos->ClkSft = ddsv;
    }
    // Ӳ���ӳ���
    pos->Tgd1 = eph->TGD1;
    pos->Tgd2 = eph->TGD2;

    return true;
}

/*
*********************************************************************
�������������źŷ���ʱ�̵�����λ�á��ٶȡ��Ӳ����
������obs     ��ǰ��Ԫ�۲�����
      GPSEph  GPS����
      BDSEph  BDS����
      xyz     ���ջ���XYZ����(m)
�������ܣ������źŷ���ʱ�̵�����λ�á��ٶȡ��Ӳ���٣�������λ�ú��ٶȣ�����߶Ƚǡ���λ�ǡ��������ӳ�
*********************************************************************
*/
void CalSatPos(epoch_t* obs, gpseph_t* GPSEph, gpseph_t* BDSEph, XYZ xyz)
{
    GPSTIME t_clock = obs->Time;
    for (int i = 0; i < obs->SatNum; i++)
    {
        GPSTIME t_AtSignalTrans;
        double dt = 0;// �������������Ӳ���ֵ�ж�
        if (obs->SatObs[i].Sys == sys_t::GPS)
        {
            gpseph_t* eph = GPSEph + obs->SatObs[i].Prn - 1;// ȡ�ñ���۲�ֵ��Ӧ�����ǵ�����
            obs->SatPVT[i].ClkOft = 0;// ��ʼ�Ӳ�Ϊ0
            do
            {   // ���������źŷ���ʱ��
                t_AtSignalTrans.Week = t_clock.Week;
                t_AtSignalTrans.SecOfWeek = t_clock.SecOfWeek - obs->SatObs[i].p[0] / CLIGHT- obs->SatPVT[i].ClkOft;
                // ���������Ӳ�
                double t_toc = diffTime(&t_AtSignalTrans, &eph->TOC);
                double st_tmp = eph->ClkBias + eph->ClkDrift * t_toc + eph->ClkDriftRate * t_toc * t_toc;
                // �����ж�ֵ
                dt = fabs(st_tmp - obs->SatPVT[i].ClkOft);
                // ���������Ӳ�
                obs->SatPVT[i].ClkOft = st_tmp;
            } while (dt > 1e-12);
            // �����Ӳ�����������
            // ��������λ�á��ٶȡ��Ӳ������
            if (CalGPSPos(obs->SatObs[i].Prn, &t_AtSignalTrans, eph, obs->SatPVT + i))
            {
                // �����źŴ���ʱ��
                double t_trans = sqrt((obs->SatPVT[i].Pos[0] - xyz.x) * (obs->SatPVT[i].Pos[0] - xyz.x) + (obs->SatPVT[i].Pos[1] - xyz.y) * (obs->SatPVT[i].Pos[1] - xyz.y) + (obs->SatPVT[i].Pos[2] - xyz.z) * (obs->SatPVT[i].Pos[2] - xyz.z)) / CLIGHT;
                // ������ת����
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

                // λ�ø���
                Matrix pos_new = Rz * pos;
                obs->SatPVT[i].Pos[0] = pos_new(0, 0);
                obs->SatPVT[i].Pos[1] = pos_new(1, 0);
                obs->SatPVT[i].Pos[2] = pos_new(2, 0);
                // �ٶȸ���
                Matrix vel_new = Rz * vel;
                obs->SatPVT[i].V[0] = vel_new(0, 0);
                obs->SatPVT[i].V[1] = vel_new(1, 0);
                obs->SatPVT[i].V[2] = vel_new(2, 0);

                // �������ǵĸ߶ȽǺͷ�λ��
                BLH blh;
                XYZToBLH(xyz, &blh, R_WGS84, F_WGS84);
                // �����վ��ƽ����ת������
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

                // ������������
                obs->SatPVT[i].TropCorr = hopfield(blh.H, obs->SatPVT[i].Elevation);
            }
        }
        else if (obs->SatObs[i].Sys == sys_t::BDS)
        {

            gpseph_t* eph = BDSEph + obs->SatObs[i].Prn - 1;// ȡ�ñ���۲�ֵ��Ӧ�����ǵ�����
            obs->SatPVT[i].ClkOft = 0;// ��ʼ�Ӳ�Ϊ0
            do
            {   // ���������źŷ���ʱ��
                t_AtSignalTrans.Week = t_clock.Week;
                t_AtSignalTrans.SecOfWeek = t_clock.SecOfWeek - obs->SatObs[i].p[0] / CLIGHT - obs->SatPVT[i].ClkOft;
                // ���������Ӳ�
                GPSTIME tb;
                tb.Week = t_AtSignalTrans.Week - 1356;
                tb.SecOfWeek = t_AtSignalTrans.SecOfWeek - 14;
                double t_toc = diffTime(&tb, &eph->TOC);
                double st_tmp = eph->ClkBias + eph->ClkDrift * t_toc + eph->ClkDriftRate * t_toc * t_toc;
                // �����ж�ֵ
                dt = fabs(st_tmp - obs->SatPVT[i].ClkOft);
                // ���������Ӳ�
                obs->SatPVT[i].ClkOft = st_tmp;
            } while (dt > 1e-12);
            // �����Ӳ�����������
            // ��������λ�á��ٶȡ��Ӳ������
            if (CalBDSPos(obs->SatObs[i].Prn, &t_AtSignalTrans, eph, obs->SatPVT + i))
            {
                // �����źŴ���ʱ��
                double t_trans = sqrt((obs->SatPVT[i].Pos[0] - xyz.x) * (obs->SatPVT[i].Pos[0] - xyz.x) + (obs->SatPVT[i].Pos[1] - xyz.y) * (obs->SatPVT[i].Pos[1] - xyz.y) + (obs->SatPVT[i].Pos[2] - xyz.z) * (obs->SatPVT[i].Pos[2] - xyz.z)) / CLIGHT;
                // ������ת����
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

                // λ�ø���
                Matrix pos_new = Rz * pos;
                obs->SatPVT[i].Pos[0] = pos_new(0, 0);
                obs->SatPVT[i].Pos[1] = pos_new(1, 0);
                obs->SatPVT[i].Pos[2] = pos_new(2, 0);
                // �ٶȸ���
                Matrix vel_new = Rz * vel;
                obs->SatPVT[i].V[0] = vel_new(0, 0);
                obs->SatPVT[i].V[1] = vel_new(1, 0);
                obs->SatPVT[i].V[2] = vel_new(2, 0);

                // �������ǵĸ߶ȽǺͷ�λ��
                BLH blh;
                XYZToBLH(xyz, &blh, R_CGS2K, F_CGS2K);
                // �����վ��ƽ����ת������
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

                // ������������
                obs->SatPVT[i].TropCorr = hopfield(blh.H, obs->SatPVT[i].Elevation);
            }
        }
    }
}