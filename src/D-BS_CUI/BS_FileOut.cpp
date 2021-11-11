/*******************************************************************************
!								"BS_FileOut.cpp"
!													2020/02/10	[Core-T]	楠崎

!*******************************************************************************/

#include "BS_FileOut.h"


int BS_FileOut::write_precision;
int BS_FileOut::iRG;
int BS_FileOut::iCY;
int BS_FileOut::iBCP;
int BS_FileOut::iBBP;

void BS_FileOut::Nothing(void) {
}


// 入力を読み込み，書き込みに必要な変数を確保しておく．
void BS_FileOut::init(const BS_FileIn&FI) {

	this->allocate(FI);

	BS_FileOut::write_precision = 15;
	this->ncn = FI.nut.size();	// ナット分割数．
	this->NT_CY.resize(this->ncn);

	this->ncs = 1;			// ひとまずのシャフト分割数．
	this->ST_CY.resize(this->ncs);

	BS_FileOut::iRG = 22;	// 剛体1つの標準出力項目数
	this->BL_ = VectorXd(BS_FileOut::iRG);
	this->NT_ = VectorXd(BS_FileOut::iRG);
	this->ST_ = VectorXd(BS_FileOut::iRG);

	BS_FileOut::iCY = 28;	// シャフト/ナット1個の出力項目数
	this->NT_CY_ = VectorXd(BS_FileOut::iCY);
	this->ST_CY_ = VectorXd(BS_FileOut::iCY);

	BS_FileOut::iBCP = 43;	// 溝-玉1組の接触出力項目数
	this->BNP_ = VectorXd(BS_FileOut::iBCP);
	this->BSP_ = VectorXd(BS_FileOut::iBCP);

	BS_FileOut::iBBP = 36;	// 玉-玉1組の接触出力項目数
	this->BBP_ = VectorXd(BS_FileOut::iBBP);

	this->np = 0;
	for (size_t i = 0; i < FI.circuit.size(); i++)
		this->np += FI.circuit[i].ball.size();

	this->CC.resize(FI.circuit.size());
	for (size_t i = 0; i < FI.circuit.size(); i++)
		this->CC[i].BL.resize(FI.circuit[i].ball.size());

	this->BNP.resize(this->np);
	this->BSP.resize(this->np);

	this->BBP.resize(this->np);

	// 時刻tの分だけ１つ多く配列長を設定
	this->BL_All = VectorXd(BS_FileOut::iRG * this->np + 1);

	this->NT_All = VectorXd(BS_FileOut::iRG + BS_FileOut::iCY * this->ncn + 1);
	this->ST_All = VectorXd(BS_FileOut::iRG + BS_FileOut::iCY * this->ncs + 1);

	this->BNP_All = VectorXd(BS_FileOut::iBCP * this->np + 1);
	this->BSP_All = VectorXd(BS_FileOut::iBCP * this->np + 1);

	this->BBP_All = VectorXd(BS_FileOut::iBBP * this->np + 1);

	this->_01_Ball_csv = FI.output._01 + ".csv";
	this->_02_Nut_csv = FI.output._02 + ".csv";
	this->_03_Shaft_csv = FI.output._03 + ".csv";
	this->_04_BallNutPair_csv = FI.output._04 + ".csv";
	this->_05_BallShaftPair_csv = FI.output._05 + ".csv";
	this->_06_BallBallPair_csv = FI.output._06 + ".csv";

	return;
}

// FileInから読み取った設定をもとに，各csvファイルのヘッダーのみを書き出すメソッド．
void BS_FileOut::write_AllHeader(void) {

	// 01_Ball.csv を書き出し．
	string ball_header = "Time[ms],";

	for (size_t i = 0; i < this->CC.size(); i++) {
		string si = to_string(i);
		for (size_t j = 0; j < this->CC[i].BL.size(); j++) {
			string si_sj = si + "_" + to_string(j);
			ball_header
				+= "xx" + si_sj + "[mm],xy" + si_sj + "[mm],xz" + si_sj + "[mm],vx" + si_sj + "[m/s],vy" + si_sj + "[m/s],vz" + si_sj + "[m/s],qx" + si_sj + "[-],qy" + si_sj + "[-],qz" + si_sj + "[-],qw" + si_sj + "[-],wx" + si_sj + "[rpm],wy" + si_sj + "[rpm],wz" + si_sj + "[rpm],axisx" + si_sj + "[-],axisy" + si_sj + "[-],axisz" + si_sj + "[-],Fx" + si_sj + "[N],Fy" + si_sj + "[N],Fz" + si_sj + "[N],Tx" + si_sj + "[Nmm],Ty" + si_sj + "[Nmm],Tz" + si_sj + "[Nmm],";
		}
	}
	ball_header += "\n";
	FileOut::write_header(this->_01_Ball_csv, ball_header);

	string nut_header = "Time[ms],xx_all[mm],xy_all[mm],xz_all[mm],vx_all[m/s],vy_all[m/s],vz_all[m/s],qx_all[-],qy_all[-],qz_all[-],qw_all[-],wx_all[rpm],wy_all[rpm],wz_all[rpm],axisx_all[-],axisy_all[-],axisz_all[-],Fx_all[N],Fy_all[N],Fz_all[N],Tx_all[Nmm],Ty_all[Nmm],Tz_all[Nmm],";
	for (int i = 0; i < this->ncn; i++) {
		string si = to_string(i);
		nut_header += "xx" + si + "[mm],xy" + si + "[mm],xz" + si + "[mm],vx" + si + "[m/s],vy" + si + "[m/s],vz" + si + "[m/s],qx" + si + "[-],qy" + si + "[-],qz" + si + "[-],qw" + si + "[-],wx" + si + "[rpm],wy" + si + "[rpm],wz" + si + "[rpm],axisx" + si + "[-],axisy" + si + "[-],axisz" + si + "[-],Fx" + si + "[N],Fy" + si + "[N],Fz" + si + "[N],Tx" + si + "[Nmm],Ty" + si + "[Nmm],Tz" + si + "[Nmm],Fsx" + si + "[N],Fsy" + si + "[N],Fsz" + si + "[N],Tsx" + si + "[Nmm],Tsy" + si + "[Nmm],Tsz" + si + "[Nmm],";
	}
	nut_header += "\n";
	FileOut::write_header(this->_02_Nut_csv, nut_header);

	string shaft_header = "Time[ms],xx_all[mm],xy_all[mm],xz_all[mm],vx_all[m/s],vy_all[m/s],vz_all[m/s],qx_all[-],qy_all[-],qz_all[-],qw_all[-],wx_all[rpm],wy_all[rpm],wz_all[rpm],axisx_all[-],axisy_all[-],axisz_all[-],Fx_all[N],Fy_all[N],Fz_all[N],Tx_all[Nmm],Ty_all[Nmm],Tz_all[Nmm],";
	for (int i = 0; i < this->ncs; i++) {
		string si = to_string(i);
		shaft_header += "xx" + si + "[mm],xy" + si + "[mm],xz" + si + "[mm],vx" + si + "[m/s],vy" + si + "[m/s],vz" + si + "[m/s],qx" + si + "[-],qy" + si + "[-],qz" + si + "[-],qw" + si + "[-],wx" + si + "[rpm],wy" + si + "[rpm],wz" + si + "[rpm],axisx" + si + "[-],axisy" + si + "[-],axisz" + si + "[-],Fx" + si + "[N],Fy" + si + "[N],Fz" + si + "[N],Tx" + si + "[Nmm],Ty" + si + "[Nmm],Tz" + si + "[Nmm],Fsx" + si + "[N],Fsy" + si + "[N],Fsz" + si + "[N],Tsx" + si + "[Nmm],Tsy" + si + "[Nmm],Tsz" + si + "[Nmm],";
	}
	shaft_header += "\n";
	FileOut::write_header(this->_03_Shaft_csv, shaft_header);

	string cylinderpair_header = "Time[ms],";
	for (int i = 0; i < this->np; i++) {
		string si = to_string(i);
		cylinderpair_header
			+= "theta" + si + "[deg],eta" + si + "[mm],zeta" + si + "[mm],";
		for (int j = 0; j < 2; j++) {
			string sij = to_string(i) + "-" + to_string(j);
			cylinderpair_header
				+= "Fnxai" + sij + "[N],Fneta" + sij + "[N],Fnzeta" + sij + "[N],"
				+ "Fsxai" + sij + "[N],Fseta" + sij + "[N],Fszeta" + sij + "[N],"
				+ "usxai" + sij + "[m/s],useta" + sij + "[m/s],uszeta" + sij + "[m/s],"
				+ "urxai" + sij + "[m/s],ureta" + sij + "[m/s],urzeta" + sij + "[m/s],"
				+ "dx" + sij
				+ "[mm],phi" + sij
				+ "[deg],a" + sij
				+ "[mm],b" + sij
				+ "[mm],h" + sij
				+ "[mm],fratio" + sij
				+ "[-],lambda" + sij
				+ "[mm],Pmax" + sij
				+ "[MPa],";
		}
	}
	cylinderpair_header += "\n";
	FileOut::write_header(this->_04_BallNutPair_csv, cylinderpair_header);
	FileOut::write_header(this->_05_BallShaftPair_csv, cylinderpair_header);
	FileOut::write_header(this->_06_BallBallPair_csv, cylinderpair_header);

	return;
}

// 各パラメータをファイル出力
void BS_FileOut::write_AllParams(double Time) {
	double t = Unit::s2ms(Time);
	int iter = 0; 
	this->BL_All[0] = t;
	for (size_t i = 0; i < this->CC.size(); i++) {
		for (size_t j = 0; j < this->CC[i].BL.size(); j++) {
			this->BL_ = BS_FileOut::form_RG(this->CC[i].BL[j]);
			for (int k = 0; k < BS_FileOut::iRG; k++) {
				int n = BS_FileOut::iRG * iter;
				this->BL_All[n + k + 1] = this->BL_[k];
			}
			iter++;
		}
	}
	this->NT_ = BS_FileOut::form_RG(this->NT);
	this->ST_ = BS_FileOut::form_RG(this->ST);

	this->NT_All[0] = t;
	this->ST_All[0] = t;
	for (int i = 0; i < BS_FileOut::iRG; i++) {
		this->NT_All[i + 1] = this->NT_[i];
		this->ST_All[i + 1] = this->ST_[i];
	}
	for (int i = 0; i < this->ncn; i++) {
		int i1 = BS_FileOut::iRG + i * BS_FileOut::iCY;
		this->NT_CY_ = BS_FileOut::form_CY(this->NT_CY[i]);
		for (int j = 0; j < BS_FileOut::iCY; j++)
			this->NT_All[i1 + j + 1] = this->NT_CY_[j];
	}
	for (int i = 0; i < this->ncs; i++) {
		int i1 = BS_FileOut::iRG + i * BS_FileOut::iCY;
		this->ST_CY_ = BS_FileOut::form_CY(this->ST_CY[i]);
		for (int j = 0; j < BS_FileOut::iCY; j++)
			this->ST_All[i1 + j + 1] = this->ST_CY_[j];
	}
	this->BNP_All[0] = t;
	this->BSP_All[0] = t;
	for (int i = 0; i < this->np; i++) {
		this->BNP_ = BS_FileOut::form_BCP(this->BNP[i]);
		this->BSP_ = BS_FileOut::form_BCP(this->BSP[i]);
		for (int j = 0; j < BS_FileOut::iBCP; j++) {
			int n = BS_FileOut::iBCP * i;
			this->BNP_All[n + j + 1] = this->BNP_[j];
			this->BSP_All[n + j + 1] = this->BSP_[j];
		}
	}
	this->BBP_All[0] = t;
	for (int i = 0; i < this->np; i++) {
		this->BBP_ = BS_FileOut::form_BBP(this->BBP[i]);
		for (int j = 0; j < BS_FileOut::iBBP; j++) {
			int n = BS_FileOut::iBBP * i;
			this->BBP_All[n + j + 1] = this->BBP_[j];
		}
	}
	// 各.csv を書き出し．
	FileOut::write_vector(this->_01_Ball_csv, this->BL_All, this->write_precision);
	FileOut::write_vector(this->_02_Nut_csv, this->NT_All, this->write_precision);
	FileOut::write_vector(this->_03_Shaft_csv, this->ST_All, this->write_precision);
	FileOut::write_vector(this->_04_BallNutPair_csv, this->BNP_All, this->write_precision);
	FileOut::write_vector(this->_05_BallShaftPair_csv, this->BSP_All, this->write_precision);
	FileOut::write_vector(this->_06_BallBallPair_csv, this->BBP_All, this->write_precision);

	return;
}

VectorXd BS_FileOut::form_RG(
	const BS_Out::Rigid & RG	// in
) {
	VectorXd RG_(BS_FileOut::iRG);
	Vector3d x_ = Unit::m2mm(Vector3d(RG.x));
	Vector3d w_ = Unit::rads2rpm(Vector3d(RG.w));
	Vector3d T_ = Unit::Nm2Nmm(Vector3d(RG.T));

	for (int i = 0; i < 3; i++) {
		RG_[0 + i] = x_[i];
		RG_[3 + i] = RG.v[i];
		RG_[10 + i] = w_[i];
		RG_[13 + i] = RG.ax[i];
		RG_[16 + i] = RG.F[i];
		RG_[19 + i] = T_[i];
	}
	for (int i = 0; i < 4; i++)
		RG_[6 + i] = RG.q[i];

	return RG_;
}

VectorXd BS_FileOut::form_CY(
	const BS_Out::Cylinder & CY	// in
) {
	VectorXd CY_(BS_FileOut::iCY);
	Vector3d x_ = Unit::m2mm(Vector3d(CY.x));
	Vector3d w_ = Unit::rads2rpm(Vector3d(CY.w));
	Vector3d T_ = Unit::Nm2Nmm(Vector3d(CY.T));
	Vector3d Ts_ = Unit::Nm2Nmm(Vector3d(CY.Ts));

	for (int i = 0; i < 3; i++) {
		CY_[0 + i] = x_[i];
		CY_[3 + i] = CY.v[i];
		CY_[10 + i] = w_[i];
		CY_[13 + i] = CY.ax[i];
		CY_[16 + i] = CY.F[i];
		CY_[19 + i] = T_[i];
		CY_[22 + i] = CY.Fs[i];
		CY_[25 + i] = Ts_[i];
	}
	for (int i = 0; i < 4; i++)
		CY_[6 + i] = CY.q[i];

	return CY_;
}

VectorXd BS_FileOut::form_BCP(
	const BS_Out::BallCylinderPair&BCP
) {
	VectorXd BCP_(BS_FileOut::iBCP);
	BCP_[0] = Unit::rad2deg(BCP.eta[0]);
	BCP_[1] = Unit::m2mm(BCP.eta[1]);
	BCP_[2] = Unit::m2mm(BCP.eta[2]);

	for (int i = 0; i < 2; i++) {
		int n = (BS_FileOut::iBCP - 3) / 2 * i;
		for (int j = 0; j < 3; j++) {
			BCP_[n + 3 + j] = BCP.GV[i].Fn[j];
			BCP_[n + 6 + j] = BCP.GV[i].Fs[j];
			BCP_[n + 9 + j] = BCP.GV[i].us[j];
			BCP_[n + 12 + j] = BCP.GV[i].ur[j];
		}
		BCP_[n + 15] = Unit::m2mm(BCP.GV[i].dx);
		BCP_[n + 16] = Unit::rad2deg(BCP.GV[i].phi);
		BCP_[n + 17] = Unit::m2mm(BCP.GV[i].a);
		BCP_[n + 18] = Unit::m2mm(BCP.GV[i].b);
		BCP_[n + 19] = Unit::m2mm(BCP.GV[i].h);
		BCP_[n + 20] = BCP.GV[i].fratio;
		BCP_[n + 21] = BCP.GV[i].lambda;
		BCP_[n + 22] = Unit::Pa2MPa(BCP.GV[i].Pmax);
	}
	return BCP_;
}

// 玉-玉の接触結果を出力．（現在はすべて0が出力されるように設定）
VectorXd BS_FileOut::form_BBP(
	const BS_Out::BallBallPair&BBP
) {
	VectorXd BBP_(BS_FileOut::iBBP);
	BBP_ = VectorXd::Zero(BS_FileOut::iBBP);
	return BBP_;
}

// 最新の数値をtempファイルから読み込む
void BS_FileOut::read_Params(
	string name,		// in : 一時ファイルの名前
	int nX,				// in : 無次元変数yの配列長
	double&t,			// out: 計算できた最後の時間
	double*Params		// out: 無次元変数yの配列
) {
	//int n = nX + 1; // tempファイルの列の長さ．時間tの分の長さを加えておく
	cout << nX << endl;
	// tempファイルは100万行は超えないと予想されるため，最後の行まで一行一行読み取る
	ifstream ifs_(name, std::ios::in);
	if (!ifs_) {
		cout << "Tempファイルがひらけません．" << endl;
	}
	string line0, line1;
	while (getline(ifs_, line0)) {
		vector<string> strvec = FileIn::split(line0, ',');
		if (strvec.size() == nX + 1) {
			line1 = line0;
		}
	}

	ifs_.close();
	vector<string> strvec = FileIn::split(line1, ',');
	cout << strvec.size() << endl;
	t = stod(strvec.at(0));

	for (int i = 0; i < nX; i++) {
		Params[i] = stod(strvec.at(i + 1));
	}

	return;
};





