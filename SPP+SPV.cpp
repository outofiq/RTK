#include "header.h"

/*
**************************************************************************
�����������㶨λ
������obs     ��ǰ��Ԫ�Ĺ۲�����
      GPSEph  GPS����
      BDSEph  BDS����
      pos     �û���λ���
             ������ʱΪ��һ��Ԫ�Ķ�λ������������к�Ϊ������Ԫ�Ķ�λ�����
����ֵ�����㶨λ�Ƿ���� true-�ɹ� false-ʧ��
�������ܣ�����SPP���㲢���������pos
**************************************************************************
*/
bool SPP(epoch_t* obs, gpseph_t* GPSEph, gpseph_t* BDSEph, epochPOS_t* pos,CFGINFO_t cfg)
{
    pos->Time = obs->Time;
    // �趨��ʼλ��
    // ��һ����Ԫ�ĳ�ʼλ������Ϊ0�����캯��������
    // ������Ԫ�ĳ�ʼλ������Ϊ��һ��Ԫ������
    XYZ X0;
    X0.x = pos->Pos[0];
    X0.y = pos->Pos[1];
    X0.z = pos->Pos[2];
    double dt0[2] = { 0,0 };// ���ջ��Ӳ0��GPS��1��BDS��
    // ������־
    bool flag = true;
    // ����������
    int calcu_num = 0;
    // �������㣬ֱ������
    do 
    {
        // �����źŷ���ʱ�̵�����λ�á��Ӳ������ת�����Ͷ������ӳ�
        CalSatPos(obs, GPSEph, BDSEph, X0);
        // ���������ǵĹ۲����ݽ������Ի�
        // �Գ�ʼλ��Ϊ�ο����Թ۲ⷽ�����Ի�����B��W����ͳ�Ʋ��붨λ�ĸ�ϵͳ������������������
        int satnum[2] = { 0,0 };// ��һ��Ԫ��ΪGPS�������������ڶ���Ԫ��ΪBDS����������
         // �ܵĿ�������������
        int sum_satnum = 0;

        //ֻ����GPSϵͳ
        if (cfg.systype == 0)
        {
            Matrix B(MAXGPSPRN, 4, 0);
            Matrix w(MAXGPSPRN, 1, 0);
            Matrix P(MAXGPSPRN, MAXGPSPRN, 0);
            for (int i = 0; i < obs->SatNum; i++)
            {
                // �۲����ݲ��������дֲ����λ�ü���ʧ�ܣ������붨λ����
                if (obs->SatObs[i].Valid == false || obs->SatPVT[i].Valid == false) continue;
                if (obs->SatObs[i].Sys != sys_t::GPS) continue;
                double �� = sqrt((obs->SatPVT[i].Pos[0] - X0.x) * (obs->SatPVT[i].Pos[0] - X0.x) +
                    (obs->SatPVT[i].Pos[1] - X0.y) * (obs->SatPVT[i].Pos[1] - X0.y) +
                    (obs->SatPVT[i].Pos[2] - X0.z) * (obs->SatPVT[i].Pos[2] - X0.z));
                B(satnum[0], 0) = (X0.x - obs->SatPVT[i].Pos[0]) / ��;
                B(satnum[0], 1) = (X0.y - obs->SatPVT[i].Pos[1]) / ��;
                B(satnum[0], 2) = (X0.z - obs->SatPVT[i].Pos[2]) / ��;
                B(satnum[0], 3) = 1;
                w(satnum[0], 0) = obs->ComObs[i].PIF - (�� + dt0[0] - CLIGHT * obs->SatPVT[i].ClkOft + obs->SatPVT[i].TropCorr);
                satnum[0] ++;
                P(satnum[0], satnum[0]) = 1;
            }
            int x_num = 4;// ��������
            // ���������������Խ��㣬���˳�����
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
            // ��λ�������ۣ�����PDOP
            pos->PDOP = sqrt(N_inv(0, 0) + N_inv(1, 1) + N_inv(2, 2));
            dt0[0] = x(3, 0);
            // ��λ�������ۣ��������λȨ�����
            Matrix v = B * x - w;
            Matrix vTv = v.transpose() * P * v;
            pos->SigmaPos = sqrt(vTv(0, 0) / (satnum[0] - x_num));
            calcu_num++;
            pos->SatNum = satnum[0];
            if (calcu_num > 15) flag = false;
        }

        //ֻ����BDSϵͳ
        if (cfg.systype == 1)
        {
            Matrix B(MAXBDSPRN, 4, 0);
            Matrix w(MAXBDSPRN, 1, 0);
            Matrix P(MAXBDSPRN, MAXBDSPRN, 0);
            for (int i = 0; i < obs->SatNum; i++)
            {
                // �۲����ݲ��������дֲ����λ�ü���ʧ�ܣ������붨λ����
                if (obs->SatObs[i].Valid == false || obs->SatPVT[i].Valid == false) continue;
                if (obs->SatObs[i].Sys != sys_t::BDS) continue;
                double �� = sqrt((obs->SatPVT[i].Pos[0] - X0.x) * (obs->SatPVT[i].Pos[0] - X0.x) +
                    (obs->SatPVT[i].Pos[1] - X0.y) * (obs->SatPVT[i].Pos[1] - X0.y) +
                    (obs->SatPVT[i].Pos[2] - X0.z) * (obs->SatPVT[i].Pos[2] - X0.z));
                B(satnum[1], 0) = (X0.x - obs->SatPVT[i].Pos[0]) / ��;
                B(satnum[1], 1) = (X0.y - obs->SatPVT[i].Pos[1]) / ��;
                B(satnum[1], 2) = (X0.z - obs->SatPVT[i].Pos[2]) / ��;
                B(satnum[1], 3) = 1;
                w(satnum[1], 0) = obs->ComObs[i].PIF - (�� + dt0[1] - CLIGHT * obs->SatPVT[i].ClkOft + obs->SatPVT[i].TropCorr + CLIGHT * (FG1_BDS * FG1_BDS * obs->SatPVT[i].Tgd1) / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS));
                satnum[1] ++;
                P(satnum[1], satnum[1]) = 1;
            }
            int x_num = 4;// ��������
            // ���������������Խ��㣬���˳�����
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
            // ��λ�������ۣ�����PDOP
            pos->PDOP = sqrt(N_inv(0, 0) + N_inv(1, 1) + N_inv(2, 2));
            dt0[1] = x(3, 0);
            // ��λ�������ۣ��������λȨ�����
            Matrix v = B * x - w;
            Matrix vTv = v.transpose() * P * v;
            pos->SigmaPos = sqrt(vTv(0, 0) / (satnum[1] - x_num));
            calcu_num++;
            pos->SatNum = satnum[1];
            if (calcu_num > 15) flag = false;
        }
    
        //GPS+BDS˫ϵͳ
        if (cfg.systype == 2)
        {
            Matrix B(MAXCHANNUM, 5, 0);
            Matrix w(MAXCHANNUM, 1, 0);
            // Ȩ��
            Matrix P(MAXCHANNUM, MAXCHANNUM, 0);

            for (int i = 0; i < obs->SatNum; i++)
            {
                sum_satnum = satnum[0] + satnum[1];
                // �۲����ݲ��������дֲ����λ�ü���ʧ�ܣ������붨λ����
                if (obs->SatObs[i].Valid == false || obs->SatPVT[i].Valid == false) continue;
                double �� = sqrt((obs->SatPVT[i].Pos[0] - X0.x) * (obs->SatPVT[i].Pos[0] - X0.x) +
                    (obs->SatPVT[i].Pos[1] - X0.y) * (obs->SatPVT[i].Pos[1] - X0.y) +
                    (obs->SatPVT[i].Pos[2] - X0.z) * (obs->SatPVT[i].Pos[2] - X0.z));

                B(sum_satnum, 0) = (X0.x - obs->SatPVT[i].Pos[0]) / ��;
                B(sum_satnum, 1) = (X0.y - obs->SatPVT[i].Pos[1]) / ��;
                B(sum_satnum, 2) = (X0.z - obs->SatPVT[i].Pos[2]) / ��;
                if (obs->SatObs[i].Sys == sys_t::GPS)
                {
                    B(sum_satnum, 3) = 1;
                    B(sum_satnum, 4) = 0;
                    w(sum_satnum, 0) = obs->ComObs[i].PIF - (�� + dt0[0] - CLIGHT * obs->SatPVT[i].ClkOft + obs->SatPVT[i].TropCorr);
                    //w(sum_satnum, 0) = obs->ComObs[i].PIF - (�� + dt0[0] - CLIGHT * obs->SatPVT[i].ClkOft ); //�����ж��������
                    satnum[0]++;// GPS������������һ
                }
                else if (obs->SatObs[i].Sys == sys_t::BDS)
                {

                    B(sum_satnum, 3) = 0;
                    B(sum_satnum, 4) = 1;
                    // ע��BDS��IF��ϵĹ۲ⷽ�̻��һ��tgdӲ���ӳ�
                    w(sum_satnum, 0) = obs->ComObs[i].PIF - (�� + dt0[1] - CLIGHT * obs->SatPVT[i].ClkOft + obs->SatPVT[i].TropCorr + CLIGHT * (FG1_BDS * FG1_BDS * obs->SatPVT[i].Tgd1) / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS));
                    //w(sum_satnum, 0) = obs->ComObs[i].PIF - (�� + dt0[1] - CLIGHT * obs->SatPVT[i].ClkOft + CLIGHT * (FG1_BDS * FG1_BDS * obs->SatPVT[i].Tgd1) / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS)); //�����ж��������
                    satnum[1]++;// BDS������������һ
                }
                P(sum_satnum, sum_satnum) = 1;// ��Ȩ��
            }

            sum_satnum = satnum[0] + satnum[1];
            // ������������㣬ֱ�ӷ��ض�λʧ��
            int x_num = 3;// ��������
            for (int i = 0; i < 2; i++)
            {
                if (satnum[i] != 0) x_num++;
            }
            // ���������������Խ��㣬���˳�����
            if (sum_satnum < x_num) return false;
            // ����������N��W
            Matrix N = B.transpose() * P * B;
            Matrix W = B.transpose() * P * w;
            Matrix x(5, 1, 0);

            // GPS\BDS˫ϵͳ��������
            if (satnum[0] != 0 && satnum[1] != 0)
            {
                // ��С������� 
                Matrix N_inv = N.inverse();
                x = N_inv * W;
                double x_norm = sqrt(x(0, 0) * x(0, 0) + x(1, 0) * x(1, 0) + x(2, 0) * x(2, 0) + x(3, 0) * x(3, 0) + x(4, 0) * x(4, 0));//sign
                if (x_norm < 1e-5)  flag = false;
                //���¶�λ���
                X0.x = X0.x + x(0, 0);
                X0.y = X0.y + x(1, 0);
                X0.z = X0.z + x(2, 0);
                // ��λ�������ۣ�����PDOP
                pos->PDOP = sqrt(N_inv(0, 0) + N_inv(1, 1) + N_inv(2, 2));
                //�����Ӳ�
                dt0[0] = x(3, 0);
                dt0[1] = x(4, 0);
            }
            // ֻ��GPS��BDS����ϵͳ������
            else
            {
                // �ع�N��W
                if (satnum[0] == 0)
                {
                    // ���GPS������Ϊ0�����4�к͵�4�о�Ϊ0������ɾ��������������Ϊ4 * 4
                    N.deletRow(3);
                    N.deleteColumn(3);
                    // ��W����ɾ����4��
                    W.deletRow(3);
                    //��B����ɾ����4��
                    B.deleteColumn(3);
                }
                else if (satnum[1] == 0)
                {
                    // BDS������Ϊ0������5�к͵�5��ɾ��
                    N.deletRow(4);
                    N.deleteColumn(4);
                    // ��W����ɾ����5��
                    W.deletRow(4);
                    //��B����ɾ����5��
                    B.deleteColumn(4);
                }
                // ��С�������
                Matrix N_inv = N.inverse();
                x = N_inv * W;

                double x_norm = sqrt(x(0, 0) * x(0, 0) + x(1, 0) * x(1, 0) + x(2, 0) * x(2, 0) + x(3, 0) * x(3, 0));
                if (x_norm < 1e-5)  flag = false;
                X0.x = X0.x + x(0, 0);
                X0.y = X0.y + x(1, 0);
                X0.z = X0.z + x(2, 0);
                // ��λ�������ۣ�����PDOP
                pos->PDOP = sqrt(N_inv(0, 0) + N_inv(1, 1) + N_inv(2, 2));

                //�����Ӳ�
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

            // ��λ�������ۣ��������λȨ�����
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
���������������
������obs     ��ǰ��Ԫ�Ĺ۲�����
      pos     �û���λ���
�������ܣ�����SPV���㲢���������pos
ע��ֻ�е�SPP�ɹ����ܽ���SPV
***************************************
*/
void SPV(epoch_t* obs, epochPOS_t* pos)
{
    Matrix B(MAXCHANNUM, 4, 0);
    Matrix w(MAXCHANNUM, 1, 0);
    int valid_satnum = 0;
    for (int i = 0; i < obs->SatNum; i++)
    {
        // �۲����ݲ��������дֲ����λ�ü���ʧ�ܣ������붨λ����
        if (obs->SatObs[i].Valid == false || obs->SatPVT[i].Valid == false) continue;
        //if (obs->SatObs[i].Sys != sys_t::GPS) continue; //ֻ��GPS
        //if (obs->SatObs[i].Sys != sys_t::BDS) continue; //ֻ��BDS
        double xsr = pos->Pos[0] - obs->SatPVT[i].Pos[0];
        double ysr = pos->Pos[1] - obs->SatPVT[i].Pos[1];
        double zsr = pos->Pos[2] - obs->SatPVT[i].Pos[2];
        double �� = sqrt(xsr * xsr + ysr * ysr + zsr * zsr);
        double ��dot = -(xsr * obs->SatPVT[i].V[0] + ysr * obs->SatPVT[i].V[1] + zsr * obs->SatPVT[i].V[2]) / ��;
        B(i, 0) = xsr / ��;
        B(i, 1) = ysr / ��;
        B(i, 2) = zsr / ��;
        B(i, 3) = 1;
        w(i, 0) = obs->SatObs[i].d[0] - (��dot - CLIGHT * obs->SatPVT[i].ClkSft);
        valid_satnum++;
    }
    // ���������������Խ��㣬���˳�����
    if (valid_satnum < 4) return;

    Matrix N = B.transpose() * B;
    Matrix W = B.transpose() * w;
    // ��С������� 
    Matrix N_inv = N.inverse();
    Matrix X = N_inv * W;
    // �洢���ٽ��
    pos->Vel[0] = X(0, 0);
    pos->Vel[1] = X(1, 0);
    pos->Vel[2] = X(2, 0);
    // �������λȨ�����
    Matrix v = B * X - w;
    Matrix vTv = v.transpose() * v;
    pos->SigmaVel = sqrt(vTv(0, 0) / (valid_satnum - 4));
}

/*
*******************************************
������������������
������pos        ������
      refpos     �ο�λ��
      denu       ��λ���
�������ܣ������������ʽת������ӡ�ڿ���̨
*******************************************
*/
void OutputResult(const epochPOS_t* pos, XYZ refpos, Matrix& denu)
{
    std::cout << pos->Time.Week << " " << std::setw(7) << std::setfill(' ') << setprecision(0) << pos->Time.SecOfWeek << " ������:" << pos->SatNum << "   ����λ��:" << std::fixed << std::setprecision(3)
        << std::setw(14) << std::setfill(' ') << pos->Pos[0] << " " << std::setw(13) << std::setfill(' ') << pos->Pos[1] << " "
        << std::setw(13) << std::setfill(' ') << pos->Pos[2] << "  PDOP:" << std::setw(6) << std::setfill(' ') << pos->PDOP << "  sigmaPos:"
        << std::setw(6) << std::setfill(' ') << pos->SigmaPos << "   �����ٶ�: " << std::setw(6) << std::setfill(' ') << pos->Vel[0] << " "
        << std::setw(6) << std::setfill(' ') << pos->Vel[1] << " " << std::setw(6) << std::setfill(' ') << pos->Vel[2] << "  sigmaV:"
        << std::setw(5) << std::setfill(' ') << pos->SigmaVel << endl;

   /* // ���㶨λ���
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
    std::cout << "  ��λ���(ENU):" << std::fixed << std::setprecision(3) << std::setw(5) << std::setfill(' ') << denu(0, 0) << " "
        << std::setw(5) << std::setfill(' ') << denu(1, 0) << " "
        << std::setw(5) << std::setfill(' ') << denu(2, 0) << std::endl << std::endl;
        */
}