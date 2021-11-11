#pragma once
#include <vector>
#include "BS_In.h"

class BS_Out {

public:
	virtual void Nothing(void)=0;	// ���ۃN���X�ł��邱�Ƃ��������߁C�Ӗ��̂Ȃ��֐����`����D

public:
	// �ȉ��̒l�͑S��SI�P�ʌn�D
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
		double f_arr;		// �w���c�ڐG�׏d[N](�X�J���[)
		double fs_[3];		// ���薀�C��[N](�a���W�n)
		double ts_[3];		// ���薀�C�g���N[Nm](�a���W�n)
		double mu_cl;		// �N�[�������C�W��[-]
		double mu_tr;		// �g���N�V�����W��[-]
		double us_[3];		// ���葬�x[m/s](�a���W�n)
		double ps_[3];		// �X���C�X�В��S[m]
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



