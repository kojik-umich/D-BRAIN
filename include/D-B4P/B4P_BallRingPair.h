#pragma once
#include <Eigen\Dense>
#include <iostream>
#include "Ball.h"
#include "B4P_Ring.h"
#include "Tribology.h"
#include "B4P_In.h"
#include "B4P_Out.h"
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::ColPivHouseholderQR;
using Eigen::AngleAxisd;
using namespace std;

// �{�[��-�����O�Ԃ̃��[�v�s�ϗʂ�ۑ����Ă����N���X�D
class B4P_BallRingPair {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	Ball    *BL;
	B4P_Ring*RG;
	double sigma;	// �����O-�{�[���̍����e���D
	double mu;		// �����O-�{�[���̃N�[�������C�W���D(typical:~0.1)
	double zeta;	// �����O-�{�[���̃_���s���O�W���D(typical:~0.2)
	double K;		// �����O-�{�[���̔���`�΂˒萔�D
	double E;		// �����O-�{�[���̓��������O���D
	double lm;		// ���j�X�J�X�����D
	double m;		// ���Z���ʁD
	double lb_alpha;
	double lb_beta;
	double lb_eta;
	double lb_k;
	double clmb_s;	// �N�[�������C�̃^���W�F���g�J�[�u�̌X��
	double fh;		// �q�X�e���V�X�����W��

	Vector3d get_ur(const Vector3d&p);
	Vector3d get_us(const Vector3d&p);



	// �O���������ݗp�ϐ��Q
	Vector3d thXZ;
	struct Slice {
		// �X���C�X�Ќ���
		double f_arr;		// �w���c�ڐG�׏d[N](�X�J���[)
		Vector3d ps;		// �X���C�X�В��S(�������W�n)
		Vector3d us;		// ���葬�x[m/s](�������W�n)
		Vector3d fs;		// ���薀�C��[N](�������W�n)
		Vector3d ts;		// ���薀�C�g���N[Nm](�������W�n)
		double mu_cl;		// �N�[�������C�W��[-]
		double mu_tr;		// �g���N�V�����W��[-]
	};


	struct Groove {
		Vector3d p;			// �ڐG�_�ʒu�i�������W�n�j
		Vector3d Fn;		// �ڐG�׏d�i�������W�n�j
		Vector3d us;
		Vector3d ur;
		Vector3d Tr;		// �]����g���N�i�������W�n�j
		Vector3d Fs;
		Vector3d Ts;
		Vector3d Fb;
		Vector3d Tb;
		Vector3d Ti;
		double a;
		double b;
		double cos_alp;
		double sin_alp;
		double Pmax;		// �ő�ʈ�[Pa]
		double lambda;		// �����_�l[-]
		double fratio;		// �����ڐG����[-]

		int msmax;			// �X���C�X������[-]
		double dx;			// �ʐڋߗ�[m]

		Vector2d p_;		// �ڐG�_�ʒu�i�a���s�f�ʍ��W�n�j
		Vector3d Fn_;		// �ڐG�׏d�i�a���s�f�ʍ��W�n�j
		Vector3d us_;
		Vector3d ur_;
		Vector3d Tr_;		// �]����g���N�i�a���s�f�ʍ��W�n�j
		Vector3d Fs_;
		Vector3d Ts_;
		Vector3d Fb_;
		Vector3d Tb_;
		Vector3d Ti_;
		Slice *SL;
		// �ȉ��̕ϐ��̓������m�ۂ̊֌W��C"��O�I��"�v�Z�r���ł��p����
		Vector3d *ps;		// �e�X���C�X�̒��S�ʒu[m]
		double *ratio_slice;	// �e�X���C�X�̐ڐG�׏d���z(�X�J���[)[-]
	}GV[2];
public:
	void link(Ball*BL, B4P_Ring*RG);
	bool how_Contact(int i, const Vector3d&bl_x,  Vector3d&er, Vector3d&eg, double&dx);
	void calc_Hertz(int i, const Vector3d & bl_x, const Vector3d & er, const Vector3d & eg, double dx, double & Rx, double & Ry, Vector3d & p, double & cos_alp, double & sin_alp, double & a, double & b, double & k);
	double calc_DynamicHertz(int i, const Vector3d & bl_x, const Vector3d & bl_v, const Vector3d & er, const Vector3d & eg, double dx, double & Rx, double & Ry, Vector3d & p, double & cos_alp, double & sin_alp, double & a, double & b);
	void calc_Sliding(int i, const Vector3d & p, double a, double Pmax, double F_norm, double fratio, Vector3d & Fs, Vector3d & Tbs, Vector3d & Tis);
	void save(bool c, int i, double cos_alp, double sin_alp, double dx, double a, double b, const Vector3d & p, const Vector3d & Fn, double Pmax, const Vector3d & us, const Vector3d & ur, const Vector3d & Tr, double lambda, double fratio, const Vector3d & Fs, const Vector3d & Ts, const Vector3d & Fb, const Vector3d & Tb, const Vector3d & Ti);
	void save_Slice(int i, int j, double farr, const Vector3d& Fs_, const Vector3d& Ts_, double mu_cl, double mu_tr, const Vector3d& us, const Vector3d& ps);
	void init_Sliceparam(int i);
	virtual void calc_force(Vector3d&Fbi, Vector3d&Nbi, Vector3d&Fib, Vector3d&Nib);
	void write(B4P_Out::BallRingPair&BXP);
	void write_slice(int ig, Matrix3d xyz2XYZ, B4P_Out::BallRingPair&GV);
	void init_Tribology(const B4P_In::Tribology & FI);
	void init_Lubrication(const B4P_In::Lubrication & FI);
	void init_Slice(int msmax);
	virtual void init(const B4P_In&FI) = 0;
	virtual double ContactAngle(double cos_alp, double sin_alp) = 0;
	void get_us_ur(const Vector3d&p, Vector3d&us, Vector3d&ur);
	void get_us_ur2(const Vector3d&p, const Vector3d&er, Vector3d&us, Vector3d&ur);
	void get_Fstf(Vector3d& Fbi, Vector3d& Fib, Vector3d& Tib);
	B4P_BallRingPair();
	~B4P_BallRingPair();

	// ���I�ɐ؂�ւ��g���C�{���W�[�N���X�D
	Tribology::Hertz *HZ;
	Tribology::RollingResistance *RR;
	Tribology::FilmThickness*FT;
	Tribology::Traction*TR;
	Tribology::Coulomb*CL;
	Tribology::DampingForce*DF;
	Tribology::Hysteresis *HY;
};

class B4P_BallOuterRingPair : public B4P_BallRingPair {
public:
	void init(const B4P_In&FI);
	double ContactAngle(double cos_alp, double sin_alp);
	void set_eta_stf(const Vector3d & eta);
};

class B4P_BallInnerRingPair : public B4P_BallRingPair {
public:
	void init(const B4P_In&FI);
	double ContactAngle(double cos_alp, double sin_alp);

};



//VectorXd make_SliceParam(void);
//double get_c(double k);
		//double*fn_slice;	// �e�X���C�X�̐ڐG�׏d(�X�J���[)[N]
		//Vector3d*Ts_slice;	// �e�X���C�X�̂��ׂ薀�C���[�����g[N]
		//Vector3d *Fs_slice;	// �e�X���C�X�̊��薀�C[N]
		//double *mutr_slice;	// �e�X���C�X�̃g���N�V�����W��[-]
		//double *mucl_slice;	// �e�X���C�X�̃g���N�V�����W��[-]
		//Vector3d *us_slice;	// �e�X���C�X�̊��葬�x[m/s]
