/*******************************************************************************
!								"B4P_FileOut.cpp"
!	出力ファイルの生成を担当するクラス．出力ファイルの読み取りもこのクラスで担当
!
!													2020/02/10	[Core-T]	楠崎
!
!*******************************************************************************/

#include "B4P_FileOut.h"

int B4P_FileOut::Z;
int B4P_FileOut::Rigid_size;
int B4P_FileOut::BRP_size;
int B4P_FileOut::BCP_size;
int B4P_FileOut::Slice_size;

void B4P_FileOut::Nothing(void) {
}

// 入力を読み込み，書き込みに必要な変数を確保しておく．
void B4P_FileOut::init(const B4P_FileIn&FI) {

	B4P_FileOut::Z = FI.ballnum;
	this->msmax = FI.msmax;

	this->allocate(B4P_FileOut::Z, _MAX_CONTACT_, this->msmax);

	// スライス出力の設定
	int n = FI.OutSlice.n;
	this->OutSlice_n = n;
	this->OutSlices = new(int[n]);
	for (int i = 0; i < n; i++) {
		this->OutSlices[i] = FI.OutSlice.list[i];
	}

	// 玉1つ/内外輪1つ当たりのパラメータ数
	B4P_FileOut::Rigid_size = 22;
	B4P_FileOut::BRP_size = 111;
	B4P_FileOut::BCP_size = 20 * _MAX_CONTACT_;
	B4P_FileOut::Slice_size = this->msmax * 15 * 2;

	this->BRP_ = VectorXd(B4P_FileOut::BRP_size);
	this->BCP_ = VectorXd(B4P_FileOut::BCP_size);
	this->Rigid_ = VectorXd(B4P_FileOut::Rigid_size);
	this->Slice_ = VectorXd(B4P_FileOut::Slice_size);

	// パラメータ数 = 玉1つのパラメータ数×玉数 + 1 (1は時間tを書き込むために確保)
	int Balls_size = B4P_FileOut::Rigid_size * B4P_FileOut::Z + 1;
	int BRPs_size = B4P_FileOut::BRP_size * B4P_FileOut::Z + 1 + B4P_FileOut::Slice_size * n;
	int BCPs_size = B4P_FileOut::BCP_size * B4P_FileOut::Z + 1;
	int Ring_size = B4P_FileOut::Rigid_size + 1;


	this->Balls_ = VectorXd(Balls_size);
	this->BRPs_ = VectorXd(BRPs_size);
	this->BCPs_ = VectorXd(BCPs_size);
	this->Ring_ = VectorXd(Ring_size);

	this->_01_Ball_csv = FI.FN.Ball + ".csv";
	this->_02_Outer_csv = FI.FN.Outer + ".csv";
	this->_03_Inner_csv = FI.FN.Inner + ".csv";
	this->_04_Cage_csv = FI.FN.Cage + ".csv";
	this->_05_BallOuterPair_csv = FI.FN.BallOuterPair + ".csv";
	this->_06_BallInnerPair_csv = FI.FN.BallInnerPair + ".csv";
	this->_07_BallCagePair_csv = FI.FN.BallCagePair + ".csv";

	this->write_precision = 10;

	return;
}

// FileInから読み取った設定をもとに，各csvファイルのヘッダーのみを書き出すメソッド．
void B4P_FileOut::write_AllHeader(void) {
	// 01_Ball.csv を書き出し．
	string ball_header = "Time[ms],";
	for (int i = 0; i < this->Z; i++) {
		string si = to_string(i);
		ball_header += "xx" + si + "[mm],xy" + si + "[mm],xz" + si + "[mm],vx" + si + "[m/s],vy" + si + "[m/s],vz" + si + "[m/s],qw" + si + "[-],qx" + si + "[-],qy" + si + "[-],qz" + si + "[-],wx" + si + "[rpm],wy" + si + "[rpm],wz" + si + "[rpm],axisx" + si + "[-],axisy" + si + "[-],axisz" + si + "[-],Fx" + si + "[N],Fy" + si + "[N],Fz" + si + "[N],Tx" + si + "[Nmm],Ty" + si + "[Nmm],Tz" + si + "[Nmm],";
	}
	ball_header += "\n";
	FileOut::write_header(this->_01_Ball_csv, ball_header);

	// 02_Outer.csv, 03_Inner.csv, 04_Cage.csv を書き出し．
	string rigid_header = "Time[ms],xx[mm],xy[mm],xz[mm],vx[m/s],vy[m/s],vz[m/s],qx[-],qy[-],qz[-],qw[-],wx[rpm],wy[rpm],wz[rpm],axisx[-],axisy[-],axisz[-],Fx[N],Fy[N],Fz[N],Tx[Nmm],Ty[Nmm],Tz[Nmm],\n";
	FileOut::write_header(this->_02_Outer_csv, rigid_header);
	FileOut::write_header(this->_03_Inner_csv, rigid_header);
	FileOut::write_header(this->_04_Cage_csv, rigid_header);

	// 05_BallOuterPair.csv, 06_BallInnerPair_csv を書き出し．
	string ringpair_header = "Time[ms],";
	for (int i = 0; i < this->Z; i++) {
		string si = to_string(i);
		ringpair_header
			+= "th" + si + "[deg],X" + si + "[mm],Z" + si + "[mm],";
		for (int j = 0; j < 2; j++) {
			string sj;
			if (j == 0)
				sj = "-";
			else
				sj = "+";
			ringpair_header
				+= "px" + si + sj + "[mm],py" + si + sj + "[mm],pz" + si + sj + "[mm],"
				+ "Fnx" + si + sj + "[N],Fny" + si + sj + "[N],Fnz" + si + sj + "[N],"
				+ "Fsx" + si + sj + "[N],Fsy" + si + sj + "[N],Fsz" + si + sj + "[N],"
				+ "Tsx" + si + sj + "[Nmm],Tsy" + si + sj + "[Nmm],Tsz" + si + sj + "[Nmm],"
				+ "Frx" + si + sj + "[N],Fry" + si + sj + "[N],Frz" + si + sj + "[N],"
				+ "Trx" + si + sj + "[Nmm],Try" + si + sj + "[Nmm],Trz" + si + sj + "[Nmm],"
				+ "usx" + si + sj + "[m/s],usy" + si + sj + "[m/s],usz" + si + sj + "[m/s],"
				+ "urx" + si + sj + "[m/s],ury" + si + sj + "[m/s],urz" + si + sj + "[m/s],"
				+ "pX" + si + sj + "[mm],pZ" + si + sj + "[mm],"
				+ "FnX" + si + sj + "[N],FnY" + si + sj + "[N],FnZ" + si + sj + "[N],"
				+ "FsX" + si + sj + "[N],FsY" + si + sj + "[N],FsZ" + si + sj + "[N],"
				+ "TsX" + si + sj + "[Nmm],TsY" + si + sj + "[Nmm],TsZ" + si + sj + "[Nmm],"
				+ "FrX" + si + sj + "[N],FrY" + si + sj + "[N],FrZ" + si + sj + "[N],"
				+ "TrX" + si + sj + "[Nmm],TrY" + si + sj + "[Nmm],TrZ" + si + sj + "[Nmm],"
				+ "usX" + si + sj + "[m/s],usY" + si + sj + "[m/s],usZ" + si + sj + "[m/s],"
				+ "urX" + si + sj + "[m/s],urY" + si + sj + "[m/s],urZ" + si + sj + "[m/s],"
				+ "lambda" + si + sj + "[-],fratio" + si + sj + "[-],alp" + si + sj + "[deg]," 
				+ "a" + si + sj + "[mm],b" + si + sj + "[mm],dx" + si + sj + "[mm],Pmax" + si + sj + "[MPa],";
		}
	}
	// 出力するスライス片のヘッダーを書き出し
	for (int i = 0; i < this->OutSlice_n; i++) {
		string bi = to_string(this->OutSlices[i]);	// 玉番号
		for (int j = 0; j < 2; j++) {
			string gj;								// 溝（+x側/-x側）
			if (j == 0) {
				gj = "-";
			}
			else {
				gj = "+";
			}
			ringpair_header += this->make_slices_header(bi, gj);
		}
	}

	ringpair_header += "\n";
	FileOut::write_header(this->_05_BallOuterPair_csv, ringpair_header);
	FileOut::write_header(this->_06_BallInnerPair_csv, ringpair_header);

	// 07_BallCagePair_csv を書き出し．
	string cagepair_header = "Time[ms],";
	for (int i = 0; i < this->Z; i++) {
		for (int j = 0; j < _MAX_CONTACT_; j++) {
			string ij = to_string(i) + to_string(j);
			cagepair_header
				+= "pattern" + ij + "[-],px" + ij + "[mm],py" + ij + "[mm],pz" + ij + "[mm],"
				+ "Fnx" + ij + "[N],Fny" + ij + "[N],Fnz" + ij + "[N],"
				+ "Fsx" + ij + "[N],Fsy" + ij + "[N],Fsz" + ij + "[N],"
				+ "px'" + ij + "[mm],py'" + ij + "[mm],pz'" + ij + "[mm],"
				+ "Fnx'" + ij + "[N],Fny'" + ij + "[N],Fnz'" + ij + "[N],"
				+ "Fsx'" + ij + "[N],Fsy'" + ij + "[N],Fsz'" + ij + "[N],"
				+ "dx" + ij + "[um],";
		}
	}
	cagepair_header += "\n";
	FileOut::write_header(this->_07_BallCagePair_csv, cagepair_header);

	return;
}

// スライス片のヘッダーを出力する．ただし，同じパラメータが隣り合うようにする
// 例えば，スライス0のus_xとスライス1のus_xが隣り合うようにする．
string B4P_FileOut::make_slices_header(string bi, string gj){
	string slice_header = "";
	string ps_x, ps_y, ps_z, us_x, us_y, us_z, fs_x, fs_y, fs_z,
		ts_x, ts_y, ts_z, fn, mucl, mutr;
	for (int k = 0; k < this->msmax; k++) {
		string sk = "s" + to_string(k);
		ps_x += "pX" + bi + gj + sk + "[mm],";
		ps_y += "pY" + bi + gj + sk + "[mm],";
		ps_z += "pZ" + bi + gj + sk + "[mm],";
		us_x += "usX" + bi + gj + sk + "[m/s],";
		us_y += "usY" + bi + gj + sk + "[m/s],";
		us_z += "usZ" + bi + gj + sk + "[m/s],";
		fs_x += "fsX" + bi + gj + sk + "[N],";
		fs_y += "fsY" + bi + gj + sk + "[N],";
		fs_z += "fsZ" + bi + gj + sk + "[N],";
		ts_x += "tsX" + bi + gj + sk + "[Nm],";
		ts_y += "tsY" + bi + gj + sk + "[Nm],";
		ts_z += "tsZ" + bi + gj + sk + "[Nm],";
		fn  += "fn" + bi + gj + sk + "[N],";
		mucl += "mucl" + bi + gj + sk + "[-],";
		mutr += "mutr" + bi + gj + sk + "[-],";
	}
	slice_header += ps_x + ps_y + ps_z + us_x + us_y + us_z + fs_x + fs_y + fs_z
		+ ts_x + ts_y + ts_z + fn + mucl + mutr;
	return slice_header;
}

// すべての数値をファイルに書き出す
void B4P_FileOut::write_AllParams(
	double Time		// in: 時間[s]
) {
	double time_ms = Unit::s2ms(Time);

	// 01_Ball.csv を書き出し．
	this->Balls_[0] = time_ms;
	for (int i = 0; i < B4P_FileOut::Z; i++) {
		this->Rigid_ = B4P_FileOut::form_RG(this->BL[i]);
		int n = B4P_FileOut::Rigid_size * i;
		for (int j = 0; j < B4P_FileOut::Rigid_size; j++)
			this->Balls_[n + j + 1] = this->Rigid_[j];
	}
	FileOut::write_vector(this->_01_Ball_csv, this->Balls_, this->write_precision);

	// 02_Outer.csv を書き出し．
	this->Rigid_ = B4P_FileOut::form_RG(this->OR);
	this->Ring_ << time_ms, this->Rigid_;
	FileOut::write_vector(this->_02_Outer_csv, this->Ring_, this->write_precision);

	// 03_Inner.csv を書き出し．
	this->Rigid_ = B4P_FileOut::form_RG(this->IR);
	this->Ring_ << time_ms, this->Rigid_;
	FileOut::write_vector(this->_03_Inner_csv, this->Ring_, this->write_precision);

	// 04_Cage.csv を書き出し．
	this->Rigid_ = B4P_FileOut::form_RG(this->CG);
	this->Ring_ << time_ms, this->Rigid_;
	FileOut::write_vector(this->_04_Cage_csv, this->Ring_, this->write_precision);

	// 05_BallOuterPair.csv を書き出し．
	this->BRPs_[0] = time_ms;
	for (int i = 0; i < B4P_FileOut::Z; i++) {
		this->BRP_ = B4P_FileOut::form_BRP(this->BOP[i]);
		int n = B4P_FileOut::BRP_size * i;
		for (int j = 0; j < B4P_FileOut::BRP_size; j++)
			this->BRPs_[n + j + 1] = this->BRP_[j];
	}
	for (int i = 0; i < this->OutSlice_n; i++) {
		int bi = this->OutSlices[i];
		this->Slice_ = B4P_FileOut::form_SL(this->BOP[bi]);
		for (int j = 0; j < B4P_FileOut::Slice_size; j++) {
			int k = B4P_FileOut::BRP_size * B4P_FileOut::Z + 1 + i * B4P_FileOut::Slice_size + j;
			this->BRPs_[k] = this->Slice_[j];
		}
	}
	FileOut::write_vector(this->_05_BallOuterPair_csv, this->BRPs_, this->write_precision);

	// 06_BallInnerPair.csv を書き出し．
	this->BRPs_[0] = time_ms;
	for (int i = 0; i < B4P_FileOut::Z; i++) {
		this->BRP_ = B4P_FileOut::form_BRP(this->BIP[i]);
		int n = B4P_FileOut::BRP_size * i;
		for (int j = 0; j < B4P_FileOut::BRP_size; j++)
			this->BRPs_[n + j + 1] = this->BRP_[j];
	}
	for (int i = 0; i < this->OutSlice_n; i++) {
		int bi = this->OutSlices[i];
		this->Slice_ = B4P_FileOut::form_SL(this->BOP[bi]);
		for (int j = 0; j < B4P_FileOut::Slice_size; j++) {
			int k = B4P_FileOut::BRP_size * B4P_FileOut::Z + 1 + i * B4P_FileOut::Slice_size + j;
			this->BRPs_[k] = this->Slice_[j];
		}
	}
	FileOut::write_vector(this->_06_BallInnerPair_csv, this->BRPs_, this->write_precision);

	// 07_BallCagePair.csv を書き出し．
	this->BCPs_[0] = time_ms;
	for (int i = 0; i < B4P_FileOut::Z; i++) {
		this->BCP_ = B4P_FileOut::form_BCP(this->BCP[i]);
		int n = B4P_FileOut::BCP_size * i;
		for (int j = 0; j < B4P_FileOut::BCP_size; j++)
			this->BCPs_[n + j + 1] = this->BCP_[j];
	}
	FileOut::write_vector(this->_07_BallCagePair_csv, this->BCPs_, this->write_precision);

	return;

}

// FileInから読み取った設定をもとに，各csvファイルのヘッダーのみを書き出すメソッド．
void B4P_FileOut::write_Header(string name) {
	FileOut::write_header(name, "");
	return;
}

// 数値をファイルに書き出す
void B4P_FileOut::write_Params(
	string name,
	double Time,
	const double*Params
) {
	int n = 13 * (this->Z + 3) + 1;
	VectorXd Params_ = VectorXd(n);

	Params_(0) = Time;

	for (int i = 1; i < n; i++)
		Params_(i) = Params[i - 1];

	FileOut::write_vector(name, Params_, 15);	// double が桁落ちしない精度が15桁なので．

	return;
};


// 最新の数値をtempファイルから読み込む
void B4P_FileOut::read_Params(
	string name,		// in : 一時ファイルの名前 
	double&t,			// out: 計算できた最後の時間
	double*Params		// out: 無次元変数yの配列
) {
	int n = 13 * (this->Z + 3);
	// tempファイルは100万行は超えないと予想されるため，最後の行まで一行一行読み取る
	ifstream ifs_(name, std::ios::in);
	if (!ifs_) {
		cout << "Tempファイルがひらけません．" << endl;
	}
	string line0, line1;
	while (getline(ifs_, line0)) {
		vector<string> strvec = FileIn::split(line0, ',');
		if (strvec.size() == n + 1) {
			line1 = line0;
		}
	}
	ifs_.close();

	vector<string> strvec = FileIn::split(line1, ',');
	t = stod(strvec.at(0));

	for (int i = 0; i < n; i++) {
		Params[i] = stod(strvec.at(i + 1));
	}

	return;
};

VectorXd B4P_FileOut::form_RG(
	const B4P_Out::Rigid & RG	// in
) {
	VectorXd RG_(B4P_FileOut::Rigid_size);

	for (int i = 0; i < 3; i++) {
		RG_[0 + i] = Unit::m2mm(RG.x[i]);
		RG_[3 + i] = RG.v[i];
		RG_[10 + i] = Unit::rads2rpm(RG.w[i]);
		RG_[13 + i] = RG.ax[i];
		RG_[16 + i] = RG.F[i];
		RG_[19 + i] = Unit::Nm2Nmm(RG.T[i]);
	}
	RG_[6] = RG.q[0];
	RG_[7] = RG.q[1];
	RG_[8] = RG.q[2];
	RG_[9] = RG.q[3];

	return RG_;
}

VectorXd B4P_FileOut::form_BRP(
	const B4P_Out::BallRingPair&BRP
) {
	VectorXd BRP_(B4P_FileOut::BRP_size);

	BRP_[0] = Unit::rad2deg(BRP.th);
	BRP_[1] = Unit::m2mm(BRP.X);
	BRP_[2] = Unit::m2mm(BRP.Z);
	
	for (int i = 0; i < 2; i++) {
		int n = (B4P_FileOut::BRP_size - 3) / 2 * i + 3;

		for (int j = 0; j < 3; j++) {
			BRP_[n + 0 + j] = Unit::m2mm(BRP.GV[i].p[j]);
			BRP_[n + 3 + j] = BRP.GV[i].Fn[j];
			BRP_[n + 6 + j] = BRP.GV[i].Fs[j];
			BRP_[n + 9 + j] = Unit::Nm2Nmm(BRP.GV[i].Ts[j]);
			BRP_[n + 12 + j] = BRP.GV[i].Fr[j];
			BRP_[n + 15 + j] = Unit::Nm2Nmm(BRP.GV[i].Tr[j]);
			BRP_[n + 18 + j] = BRP.GV[i].us[j];
			BRP_[n + 21 + j] = BRP.GV[i].ur[j];
		}
		for (int j = 0; j < 2; j++)
			BRP_[n + 24 + j] = Unit::m2mm(BRP.GV[i].p_[j]);
		for (int j = 0; j < 3; j++) {
			BRP_[n + 26 + j] = BRP.GV[i].Fn_[j];
			BRP_[n + 29 + j] = BRP.GV[i].Fs_[j];
			BRP_[n + 32 + j] = Unit::Nm2Nmm(BRP.GV[i].Ts_[j]);
			BRP_[n + 35 + j] = BRP.GV[i].Fr_[j];
			BRP_[n + 38 + j] = Unit::Nm2Nmm(BRP.GV[i].Tr_[j]);
			BRP_[n + 41 + j] = BRP.GV[i].us_[j];
			BRP_[n + 44 + j] = BRP.GV[i].ur_[j];
		}
		BRP_[n + 47] = BRP.GV[i].lambda;
		BRP_[n + 48] = BRP.GV[i].fratio;
		BRP_[n + 49] = Unit::rad2deg(BRP.GV[i].phi);
		BRP_[n + 50] = Unit::m2mm(BRP.GV[i].a);
		BRP_[n + 51] = Unit::m2mm(BRP.GV[i].b);
		BRP_[n + 52] = Unit::m2mm(BRP.GV[i].dx);
		BRP_[n + 53] = Unit::Pa2MPa(BRP.GV[i].Pmax);
	}
	return BRP_;
}

VectorXd B4P_FileOut::form_BCP(
	const B4P_Out::BallCagePair&BCP
) {
	VectorXd BCP_(B4P_FileOut::BCP_size);

	for (int i = 0; i < _MAX_CONTACT_; i++) {
		int n = B4P_FileOut::BCP_size / _MAX_CONTACT_ * i;
		BCP_[n] = BCP.CP[i].ptt;
		for (int j = 0; j < 3; j++) {
			BCP_[n + 1 + j] = Unit::m2mm(BCP.CP[i].p[j]);
			BCP_[n + 4 + j] = BCP.CP[i].Fn[j];
			BCP_[n + 7 + j] = BCP.CP[i].Fs[j];
			BCP_[n + 10 + j] = Unit::m2mm(BCP.CP[i].p_[j]);
			BCP_[n + 13 + j] = BCP.CP[i].Fn_[j];
			BCP_[n + 16 + j] = BCP.CP[i].Fs_[j];
		}
		BCP_[n + 19] = Unit::m2um(BCP.CP[i].dx);
	}
	return BCP_;
}

// 各スライス片の結果を出力(スライス毎ではなく，パラメータ毎に出力)
VectorXd B4P_FileOut::form_SL(
	const B4P_Out::BallRingPair&BRP
) {
	VectorXd SL_(B4P_FileOut::Slice_size);
	for(int i = 0; i < 2; i++){
		int n = B4P_FileOut::Slice_size / 2 * i;
		for (int j = 0; j < this->msmax; j++) {
			for (int k = 0; k < 3; k++) {
				SL_[n + (0 + k) * this->msmax + j] = Unit::m2mm(BRP.GV[i].SL[j].ps_[k]);
				SL_[n + (3 + k) * this->msmax + j] = BRP.GV[i].SL[j].us_[k];
				SL_[n + (6 + k) * this->msmax + j] = BRP.GV[i].SL[j].fs_[k];
				SL_[n + (9 + k) * this->msmax + j] = BRP.GV[i].SL[j].ts_[k];
			}
			SL_[n + 12 * this->msmax + j] = BRP.GV[i].SL[j].f_arr;
			SL_[n + 13 * this->msmax + j] = BRP.GV[i].SL[j].mu_cl;
			SL_[n + 14 * this->msmax + j] = BRP.GV[i].SL[j].mu_tr;
		}

	}

	return SL_;
}