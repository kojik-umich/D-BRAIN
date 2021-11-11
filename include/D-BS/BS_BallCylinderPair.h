#pragma once
#include <Eigen\Dense>
#include "Ball.h"
#include "Tribology.h"
#include "BS_In.h"
#include "BS_Out.h"
#include "BS_Cylinder.h"
#include "BS_Nut.h"
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Matrix3d;


class BS_BallCylinderPair {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:
	Ball       *BL;		// �o�^����Ă���ʁD
	BS_Cylinder*CY;		// �o�^����Ă��闆���~�����ށD
	int         iSP;		// ��ԍ��D
	double       E;		// ���������O�� [Pa]
	double       m;		// ���Z���� [kg]

	// �ڐG�v�Z�Ɏg���g���C�{���W�[��
	struct Tribologies {
		Tribology::Hertz			*HZ;	// �w���c�ڐG��
		Tribology::RollingResistance*RR;	// �]����S����R��
		Tribology::FilmThickness	*FT;	// ����������
		Tribology::Traction			*TR;	// �����g���N�V�����͂̎�
		Tribology::Coulomb			*CL;	// �N�[�������C��
		double						cs;		// �N�[�������C�̗����オ��X���D
		Tribology::Hysteresis		*HY;	// �q�X�e���V�X����
		double						fh;		// �q�X�e���V�X�W��
	}TB;

	// �������
	struct Lubrication {
		double eta;		// �S�x [Pa*s]
		double beta;	// ���x�S�x�W�� [K^-1]
		double k;		// ���M�`���� [W/(K*m)]
		double alpha;	// ���͔S�x�W�� [Pa^-1]
		double lm;		// ���j�X�J�X���� [m]
	}LB;

	// �e�a�̕����l
	struct Groove {
		double	sigma;		// �\�ʑe�� [m]
		int		n;			// �X���C�X��
		double  *r;			// �ڐG���͂̕��z�� [-]
		Vector2d*xy;		// �X���C�X���W�� [-] xy�͒P�ʉ~�̓����̒l�D0����:a�����D1����:b�����D
		Vector3d*ps;		// �X���C�X���W [m]�i�{���̓��[�J���ϐ��Ƃ��Ĉ����������C����new����ƃR�X�g�������邽�߃����o�Ƃ��Ďg�p�j
		double	mu;			// ���C�W�� [-]
		double	zeta;		// ������ [-]
	} GV[2];
	
	// �o�͕ۑ��p�D
	struct Save {
		Vector3d eta;			// �����ʑ��p[rad], eta���W[m], zeta���W[m] �̔z��
		// �e�a�̈ꎞ�ϐ�
		struct Groove {
			Vector3d Fn;		// �����׏d [N]�i�������W�n�j
			Vector3d Fs;		// ���薀�C�i�X���C�X�̑��a�j [N]�i�������W�n�j
			Vector3d us;		// ���葬�x [m/s]�i�������W�n�j
			Vector3d ur;		// �]���葬�x [m/s]�i�������W�n�j
			double dx;			// �e���ڋߗ� [m]
			double phi;			// �ڐG�p [rad]
			double a;			// �ڐG�ȉ~���a [m]
			double b;			// �ڐG�ȉ~�Z�a [m]
			double h;			// �������� [m]
			double fratio;		// �����׏d�x�� [-]
			double lambda;		// �����_�l [-]
			Vector3d p;			// �ڐG�ȉ~�����_ [m]�i�������W�n�j
			double F;			// �S�Ă̗͂̑��a [N]
			Vector2d Fbc;		// �����׏d [N]�i�a���s�f�ʍ��W�n�j
			double Pmax;		// �ő�ʈ�[Pa]
			struct Slice {
				double f_arr;		// �w���c�ڐG�׏d[N](�X�J���[)
				Vector3d ps;		// �X���C�X�В��S(�������W�n)
				Vector3d us;		// ���葬�x[m/s](�������W�n)
				Vector3d fs;		// ���薀�C��[N](�������W�n)
				Vector3d ts;		// ���薀�C�g���N[Nm](�������W�n)
				double mu_cl;		// �N�[�������C�W��[-]
				double mu_tr;		// �g���N�V�����W��[-]
			} *SL;					// �e�X���C�X�Ђ̌��ʁi�|�C���^�z��C���� msmax�j
		} GV[2];
	} SV;

public:
	void link(Ball * BL, BS_Cylinder * CY, int is);

	bool how_Contact(int i, const Vector2d & bl_x, Vector2d & e, double & dx);
	double calc_Hertz(int ig, const Vector2d & bl_eta, const Vector2d & e, double dx, double & Rx, double & Ry, Vector2d & p, double & cos_alp, double & sin_alp, double & a, double & b, double & k, Vector2d & rho);
	double calc_DynamicHertz(int ig, const Vector2d & bl_eta, const Vector2d & bl_etav, const Vector2d & e, double dx, double & Rx, double & Ry, Vector2d & p, double & cos_alp, double & sin_alp, double & a, double & b, Vector2d & rho);
	void get_F0(Vector2d & Fbc, Vector3d & Fcb, Vector3d & Tcb);
	void save_F0(int i, const Vector3d & p, double F, const Vector3d & eta, double a, const Vector2d & Fbc);
	void get_F1(bool v, Vector2d & Fbc, Vector3d & Fcb, Vector3d & Tcb);
	void get_F2(Vector3d & vF, Vector3d & vT);
	void get_FT(Vector3d & Fbc, Vector3d & Tbc, Vector3d & Tcb, Vector3d & Fs, Vector3d & Ts);
	Vector3d get_ur(const Vector3d & p);
	Vector3d get_us(const Vector3d & p);
	Vector3d get_eta(void);
	void calc_Sliding(int ig, double th, const Vector3d & p, const Vector3d & xai, double a, double b, double Pmax, double ur_norm, double F_norm, double fratio, const Vector2d & rho, Vector3d & Fs, Vector3d & Ts, Vector3d & Tcs);
	void init(const BS_In::BallCylinderPair & BCP, const BS_In::Tribology & tribology, const BS_In::Oil & oil);
	void init_Tribology(const BS_In::Tribology & tribology);
	void init_Oil(const BS_In::Oil & oil);
	void save(int ig, const Vector3d & eta, const Vector3d & Fn, const Vector3d & Fs, const Vector3d & us, const Vector3d & ur, double dx, double cos_alp, double sin_alp, double a, double b, double h, const Vector3d & p, double lambda, double fratio, double Pmax);
	Vector3d get_e(void);
	void save(BS_Out::BallCylinderPair & OUT);
	void save_Slice(int i, int j, double farr, const Vector3d& Fs_, const Vector3d& Ts_, double mu_cl, double mu_tr, const Vector3d& us, const Vector3d& ps);
	void init_Sliceparam(int i);
	void write_slice(int ig, Matrix3d xyz2XYZ, BS_Out::BallCylinderPair&BRP);
	BS_BallCylinderPair();
	~BS_BallCylinderPair();
};


class BS_BallNutPair : public BS_BallCylinderPair {

public:
	double      th0;	// �i��ɐÉ�͂ŗ��p����j�i�b�g���猩���ʂ̈ʑ��p�D
	void set_eta0(const Vector2d&eta);
	Vector2d get_etavector0(const Vector3d & x);
};

class BS_BallShaftPair : public BS_BallCylinderPair {
public:
	double      thp;	// �i��ɐÉ�͂ŗ��p����j�i�b�g���猩���ʂ̈ʑ��p�D

	void set_etap(const Vector2d&eta);
	Vector2d get_etavectorp(const Vector3d & x);
};


//public:
//	enum Contact {
//		NC,		// �ڂ��Ă��Ȃ���ԁDNotContact�D
//		IVR,	// ��e�����̏�����ԁD���ݎ������Ă��Ȃ��D
//		EHL		// �e�����̏�����ԁD
//	};
