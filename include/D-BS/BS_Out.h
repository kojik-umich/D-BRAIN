#pragma once
#include <vector>
#include "BS_In.h"

class BS_Out {

public:
	virtual void Nothing(void)=0;	// 抽象クラスであることを示すため，意味のない関数を定義する．

public:
	// 以下の値は全てSI単位系．
	struct Rigid {
	public:
		double x[3];
		double q[4];
		double v[3];
		double w[3];
		double ax[3];
		double F[3];
		double T[3];
	} NT, ST;

	struct Circuit {
		std::vector<Rigid> BL;
	};
	std::vector<Circuit> CC;

	struct Cylinder {
	public:
		double x[3];
		double q[4];
		double v[3];
		double w[3];
		double ax[3];
		double F[3];
		double T[3];
		double Fs[3];
		double Ts[3];
	};
	std::vector<Cylinder> NT_CY;
	std::vector<Cylinder> ST_CY;

	struct Slice {
		double f_arr;		// ヘルツ接触荷重[N](スカラー)
		double fs_[3];		// 滑り摩擦力[N](溝座標系)
		double ts_[3];		// 滑り摩擦トルク[Nm](溝座標系)
		double mu_cl;		// クーロン摩擦係数[-]
		double mu_tr;		// トラクション係数[-]
		double us_[3];		// 滑り速度[m/s](溝座標系)
		double ps_[3];		// スライス片中心[m]
	};

	struct BallCylinderPair {
	public:
		double eta[3];

		struct Groove {
		public:
			double Fn[3];
			double Fs[3];
			double us[3];
			double ur[3];
			double dx;
			double phi;
			double a;
			double b;
			double h;
			double fratio;
			double lambda;
			double Pmax;
			Slice *SL;
		}GV[2];
	};
	std::vector<BallCylinderPair> BNP;
	std::vector<BallCylinderPair> BSP;

	struct BallBallPair {
	public:
		double p[3];
		double Fn[3];
		double Fs[3];
		double us[3];
		double ur[3];
		double dx[3];
		double a[3];
		double h[3];
	};
	std::vector<BallBallPair> BBP;

	void allocate(const BS_In & IN) {

		size_t nn = IN.nut.size();
		this->NT_CY.resize(nn);

		size_t ns = IN.shaft.size();
		this->ST_CY.resize(ns);

		size_t ballnum = 0;

		size_t nc = IN.circuit.size();
		this->CC.resize(nc);
		for (size_t i = 0; i < nc; i++) {
			size_t nb = IN.circuit[i].ball.size();
			this->CC[i].BL.resize(nb);
			ballnum += nb;
		}

		this->BNP.resize(ballnum);
		this->BSP.resize(ballnum);
		this->BBP.resize(ballnum);
		int msmax = 21;
		for (size_t i = 0; i < ballnum; i++) {
			for (int j = 0; j < 2; j++) {
				this->BNP[i].GV[j].SL = new Slice[msmax];
				this->BSP[i].GV[j].SL = new Slice[msmax];
			}
		}

		return;
	}
};



