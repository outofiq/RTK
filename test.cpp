#include "header.h"

unsigned char BUFF[MAXRAWLEN];                   //用于解码的数据
unsigned char buff[MAXBUFF];                     //使用网口时暂存数据的数组
gpseph_t GPSEph[MAXGPSPRN] = { gpseph_t() };     //GPS星历
gpseph_t BDSEph[MAXBDSPRN] = { gpseph_t() };     //BDS星历
epoch_t obs = epoch_t();                         //储存待解算历元观测值的结构体
epochPOS_t pos = epochPOS_t();                   //储存解算结果的结构体
Matrix dENU(3, 1, 0);                            //定位误差
XYZ filepos = XYZ();                             //文件的参考坐标
XYZ socketspos = XYZ();                          //网口SPP的参考坐标
CFGINFO_t cfg = CFGINFO_t();                     //配置参数结构体
rtkData_t rtkdata = rtkData_t();                 //rtk基准站和流动站数据
RTKResult_t result = RTKResult_t();          //储存浮点解的结构体变量

int main()
{
	//int Epoch_flag = 0;                          //解算第多少个历元
	int ReadLen = 0;                             //已经读取的数据字节数
	int judge = 0;                               //判断是否进行spp解算
	int length = 0;                              //进行读取操作时读取的字节数

	double time;

	if (!ReadCFG("CONFIG.txt", &cfg))
	{
		cerr << "读取配置文件失败！" << endl;
		return 0;
	}

	//ofstream resultfile("E:\\卫导算法\\2\\lamda_result.txt");
	ofstream resultfile("lamda_result.txt");
	if (!resultfile.is_open())
	{
		cerr << "结果文件打开失败！" << endl;
	}

	//ofstream FixResultFile("E:\\卫导算法\\2\\fixed_result.txt");
	ofstream FixResultFile("fixed_result.txt");
	if (!FixResultFile.is_open())
	{
		cerr << "结果文件打开失败" << endl;
	}

	if (cfg.datatype == 0)
	{
		//打开基准站数据文件
		ifstream FBas(cfg.BasInfo, ios::in | ios::binary);
		if (!FBas.is_open())
		{
			cerr << "打开基站文件失败！" << endl;
			return 0;
		}
		//打开流动站数据文件
		ifstream FRov(cfg.RovInfo, ios::in | ios::binary);
		if (!FRov.is_open())
		{
			cerr << "打开流动站文件失败！" << endl;
			return 0;
		}

		while (1)
		{
			//进行时间同步
			int syncflag;
			syncflag = TimeSyncPP(FBas, FRov, &rtkdata);
			if (syncflag == -1) break; // 文件结束，跳出循环
			else if (syncflag == 0) continue; // 同步失败，继续循环
			else;  //同步成功，继续解算

			DetectError(&rtkdata.basData);
			DetectError(&rtkdata.rovData);
			if (SPP(&rtkdata.basData, rtkdata.GPSEph, rtkdata.BDSEph, &rtkdata.BasPres, cfg) && SPP(&rtkdata.rovData, rtkdata.GPSEph, rtkdata.BDSEph, &rtkdata.RovPres, cfg))
			{
				SPV(&rtkdata.basData, &rtkdata.BasPres);
				SPV(&rtkdata.rovData, &rtkdata.RovPres);
				//OutputResult(&rtkdata.RovPres, filepos, dENU);
				SDEpochObs(&rtkdata.rovData, &rtkdata.basData, &rtkdata.SDObs);
				SelectRefSat(&rtkdata.rovData, &rtkdata.SDObs, &rtkdata.DDObs);
				DDEpochObs(&rtkdata.SDObs, &rtkdata.DDObs);
				if (LSRTKFloat(&rtkdata, cfg, &result))
				{
					GPSTIME time = result.Time;
					double weeks = time.SecOfWeek;
					int n = result.n;
					int m = result.m;
					double* fa = result.fa;
					double* Qa = result.Qa;
					double* F = result.F;
					double* s = result.s;

					lambda(n, 2, fa, Qa, F, s);

					double Ratio = s[1] / s[0];
					result.Ratio = Ratio;
					if (result.Ratio > 2.0) result.fixed = 1;
					else result.fixed == 0;

					for (int i = 0; i < n; i++)
					{
						result.a_fix(i, 0) = F[i];
					}

					resultfile << fixed << setprecision(4) << result.RovPos[0] << " " << result.RovPos[1] << " " << result.RovPos[2] << " ";
					resultfile << fixed << setprecision(4) << result.sigma << " " << result.mxyz[0] << " " << result.mxyz[1] << " " << result.mxyz[2] << " " << result.ms << endl;
					resultfile << std::fixed << std::setprecision(6) << weeks << "   ";
					for (int i = 0; i < n; i++) {
						resultfile << std::fixed << std::setw(15) << std::setprecision(5) << fa[i] << "	";
					}
					resultfile << "\n";

					resultfile << std::fixed << std::setprecision(6) << weeks << "   ";
					for (int i = 0; i < n; i++) {
						resultfile << std::fixed << std::setw(15) << std::setprecision(5) << F[i] << "  ";
					}
					resultfile << "\n";

					resultfile << std::fixed << std::setprecision(6) << weeks << "   ";
					for (int i = 0; i < n; i++) {
						resultfile << std::fixed << std::setw(15) << std::setprecision(5) << F[i + n] << "  ";
					}
					resultfile << "\n";

					resultfile << std::fixed << std::setprecision(6) << weeks << "   "
						<< "Sqnorm[0] = " << std::fixed << std::setw(10) << std::setprecision(5) << s[0]
						<< ", Sqnorm[1] = " << std::fixed << std::setw(10) << std::setprecision(5) << s[1]
						<< ", Ratio = " << std::fixed << std::setw(4) << std::setprecision(2) << (s[1] / s[0])
						<< "\n\n\n";

					resultfile.flush();
					free(result.fa); free(result.Qa); free(result.F); free(result.s);


					LSRTKFix(&rtkdata, &result);
					FixOutPut(FixResultFile, &result, cfg);
				}
			}
		}
		FBas.close();
		FRov.close();
		resultfile.close();
		FixResultFile.close();
		return 0;
	}


	if (cfg.datatype == 1)
	{
		SOCKET BasSock, RovSock;
		int LenRead, LenDecode = 0;
		vector<epoch_t> rovOBS, basOBS;

		//打开基准站网口
		if (OpenSocket(BasSock, cfg.BasIP, cfg.BasPort) == false)
		{
			printf("This ip & port was not opened.\n");
			return 0;
		}

		//打开流动站网口
		if (OpenSocket(RovSock, cfg.RovIP, cfg.RovPort) == false)
		{
			printf("This ip & port was not opened.\n");
			return 0;
		}

		//打开保存原始二进制数据的文件
		ofstream BasData("BasData.txt", ios::out | ios::binary);
		if (!BasData.is_open())
		{
			cerr << "open error!";
		}

		ofstream RovData("RovData.txt", ios::out | ios::binary);
		if (!RovData.is_open())
		{
			cerr << "open error!";
		}

		while (1)
		{
			//进行时间同步
			int SyncFlag;
			SyncFlag = TimeSyncRT(BasSock, RovSock, &rtkdata, rovOBS, basOBS, BasData,RovData);
			if (SyncFlag == 0) continue; //同步失败
			else;

			DetectError(&rtkdata.basData);
			DetectError(&rtkdata.rovData);
			if (SPP(&rtkdata.basData, rtkdata.GPSEph, rtkdata.BDSEph, &rtkdata.BasPres, cfg) && SPP(&rtkdata.rovData, rtkdata.GPSEph, rtkdata.BDSEph, &rtkdata.RovPres, cfg))
			{
				SPV(&rtkdata.basData, &rtkdata.BasPres);
				SPV(&rtkdata.rovData, &rtkdata.RovPres);
				//OutputResult(&rtkdata.RovPres, filepos, dENU);
				SDEpochObs(&rtkdata.rovData, &rtkdata.basData, &rtkdata.SDObs);
				SelectRefSat(&rtkdata.rovData, &rtkdata.SDObs, &rtkdata.DDObs);
				DDEpochObs(&rtkdata.SDObs, &rtkdata.DDObs);
				if (LSRTKFloat(&rtkdata, cfg, &result))
				{
					GPSTIME time = result.Time;
					double weeks = time.SecOfWeek;
					int n = result.n;
					int m = result.m;
					double* fa = result.fa;
					double* Qa = result.Qa;
					double* F = result.F;
					double* s = result.s;

					lambda(n, 2, fa, Qa, F, s);

					double Ratio = s[1] / s[0];
					result.Ratio = Ratio;
					if (result.Ratio > 2.0) result.fixed = 1;
					else result.fixed = 0;

					for (int i = 0; i < n; i++)
					{
						result.a_fix(i, 0) = F[i];
					}

					resultfile << fixed << setprecision(4) << result.RovPos[0] << " " << result.RovPos[1] << " " << result.RovPos[2] << " ";
					resultfile << fixed << setprecision(4) << result.sigma << " " << result.mxyz[0] << " " << result.mxyz[1] << " " << result.mxyz[2] << " " << result.ms << endl;
					resultfile << std::fixed << std::setprecision(6) << weeks << "   ";
					for (int i = 0; i < n; i++) {
						resultfile << std::fixed << std::setw(15) << std::setprecision(5) << fa[i] << "	";
					}
					resultfile << "\n";

					resultfile << std::fixed << std::setprecision(6) << weeks << "   ";
					for (int i = 0; i < n; i++) {
						resultfile << std::fixed << std::setw(15) << std::setprecision(5) << F[i] << "  ";
					}
					resultfile << "\n";

					resultfile << std::fixed << std::setprecision(6) << weeks << "   ";
					for (int i = 0; i < n; i++) {
						resultfile << std::fixed << std::setw(15) << std::setprecision(5) << F[i + n] << "  ";
					}
					resultfile << "\n";

					resultfile << std::fixed << std::setprecision(6) << weeks << "   "
						<< "Sqnorm[0] = " << std::fixed << std::setw(10) << std::setprecision(5) << s[0]
						<< ", Sqnorm[1] = " << std::fixed << std::setw(10) << std::setprecision(5) << s[1]
						<< ", Ratio = " << std::fixed << std::setw(4) << std::setprecision(2) << (s[1] / s[0])
						<< "\n\n\n";

					resultfile.flush();
					free(result.fa); free(result.Qa); free(result.F); free(result.s);

				}
				LSRTKFix(&rtkdata, &result);
				FixOutPut(FixResultFile, &result, cfg);
			}
		}
		CloseSocket(BasSock);
		CloseSocket(RovSock);
		resultfile.close();
		FixResultFile.close();
		BasData.close();
		RovData.close();
		return 0;
	}
}