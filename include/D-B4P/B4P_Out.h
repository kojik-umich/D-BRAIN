#pragma once

class B4P_Out {
public:
	virtual void Nothing(void) = 0;	// 抽象クラスであることを示すため，意味のない関数を定義する．

public:
	// 以下の値は全てSI単位系．
	struct Rigid {
		double x[3];
		double q[4];	// q[0]:x成分，q[1]:y成分，q[2]:z成分，q[3]:w成分
		double v[3];
		double w[3];
		double ax[3];
		double F[3];
		double T[3];
	} *BL, OR, IR, CG;

	struct Slice {
		double f_arr;		// ヘルツ接触荷重[N](スカラー)
		double fs_[3];		// 滑り摩擦力[N](溝座標系)
		double ts_[3];		// 滑り摩擦トルク[Nm](溝座標系)
		double mu_cl;		// クーロン摩擦係数[-]
		double mu_tr;		// トラクション係数[-]
		double us_[3];		// 滑り速度[m/s](溝座標系)
		double ps_[3];		// スライス片中心[m]
	};

	struct Groove {
		double p[3];
		double Fn[3];
		double Fs[3];
		double Ts[3];
		double Fr[3];
		double Tr[3];
		double us[3];
		double ur[3];
		double fratio;
		double lambda;
		double phi;
		double a;
		double b;
		double dx;
		double Pmax;

		double p_[3];
		double Fn_[3];
		double Fs_[3];
		double Ts_[3];
		double Fr_[3];
		double Tr_[3];
		double us_[3];
		double ur_[3];
		Slice *SL;
	};


	struct BallRingPair {
		double th;
		double X;
		double Z;
		Groove GV[2];
	}*BOP, *BIP;
 

	struct BallCagePair {
		struct ContactPoint {
			double p[3];
			double Fn[3];
			double Fs[3];
			double p_[3];
			double Fn_[3];
			double Fs_[3];
			double dx;
			int ptt;
		} *CP;
	} *BCP;

	B4P_Out() {
		this->BL = NULL;
		this->BOP = NULL;
		this->BIP = NULL;
		this->BCP = NULL;
		return;
	};
	// メンバ変数の配列の動的確保
	void allocate(
		int Z,		// in: 玉数
		int m,		// in: 玉ー保持器接触点数
		int msmax	// in: スライス片個数
	) {
		this->BL = new Rigid[Z];
		this->BOP = new BallRingPair[Z];
		this->BIP = new BallRingPair[Z];
		this->BCP = new BallCagePair[Z];

		for (int i = 0; i < Z; i++) {
			this->BCP[i].CP = new BallCagePair::ContactPoint[m];
			for (int j = 0; j < 2; j++) {
				this->BOP[i].GV[j].SL = new Slice[msmax];
				this->BIP[i].GV[j].SL = new Slice[msmax];
			}
		}
		return;
	};


	~B4P_Out() {
		if (this->BL != NULL)
			delete[] this->BL;
		if (this->BOP != NULL)
			delete[] this->BOP;
		if (this->BIP != NULL)
			delete[] this->BIP;
		if (this->BCP != NULL)
			delete[] this->BCP;
		return;
	};
};



