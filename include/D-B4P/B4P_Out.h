#pragma once

class B4P_Out {
public:
	virtual void Nothing(void) = 0;	// ���ۃN���X�ł��邱�Ƃ��������߁C�Ӗ��̂Ȃ��֐����`����D

public:
	// �ȉ��̒l�͑S��SI�P�ʌn�D
	struct Rigid {
		double x[3];
		double q[4];	// q[0]:x�����Cq[1]:y�����Cq[2]:z�����Cq[3]:w����
		double v[3];
		double w[3];
		double ax[3];
		double F[3];
		double T[3];
	} *BL, OR, IR, CG;

	struct Slice {
		double f_arr;		// �w���c�ڐG�׏d[N](�X�J���[)
		double fs_[3];		// ���薀�C��[N](�a���W�n)
		double ts_[3];		// ���薀�C�g���N[Nm](�a���W�n)
		double mu_cl;		// �N�[�������C�W��[-]
		double mu_tr;		// �g���N�V�����W��[-]
		double us_[3];		// ���葬�x[m/s](�a���W�n)
		double ps_[3];		// �X���C�X�В��S[m]
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
	// �����o�ϐ��̔z��̓��I�m��
	void allocate(
		int Z,		// in: �ʐ�
		int m,		// in: �ʁ[�ێ���ڐG�_��
		int msmax	// in: �X���C�X�Ќ�
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



