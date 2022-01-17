#pragma once
#include <stdio.h>
#include <Eigen\Dense>
#include "BS_BallScrew.h"
#include "intel_ode.h"
#include "mkl_rci.h"
#include "Rigid.h"
#include "BS_In.h"
#include "BS_FileIn.h"

class BS_Calculator {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
public:
	//static BS_BallScrew BS;
	static std::shared_ptr<BS_BallScrew> BS;

	// �É�͂̂��߂ɕK�v�ȃO���[�o���ϐ��p�����^�Q�DMKL�ł͊֐��̌`�����܂��Ă��邽�߁C���̂悤�Ȍ`�Œl���󂯓n���D
	static struct Stt {
		struct Set {
			int n;				// ���͔z��T�C�Y�D
			int m;				// �o�͔z��T�C�Y�D
			double eps[6];		// �\���o�̐ݒ�i����6�̔z��j�D�e�v�f�̈Ӗ��̓}�j���A���Q�ƁD
			int iter1;			// �ő�J��Ԃ����D
			int iter2;			// �ő厎�s�񐔁D
			double rs;			// �����X�e�b�v�����D
			double jac_eps;		// ���R�r�A���v�Z�̐��x�D
		} set[3];
		double v0;		// [m/s]:	�V���t�g���i���x�i//axis�j
		double w0;		// [rad/s]:	�V���t�g��]���x�i//axis�j
		double wn;		// [rad/s]:	�i�b�g��]���x�i//axis�j
		int i1;			// �X�e�b�v1�p�D�ʔԍ����i�[�D
	} stt;

	static struct Dyn {
		struct Set {
			double t_end;		// �v�Z�I������[s]
			double t_step;		// �T���v�����O���ԃX�e�b�v�T�C�Y�i���ŏ����ԃX�e�b�v�T�C�Y�j
			int nX;				// �ʁE���O�ցE�ێ���e�v�f�̏�ԗʁi���W�E���x�E�����x�j
			double h;			// �ŏ��X�e�b�v�T�C�Y
			double hmin;		// ���e�덷
			double ep;			// �������l
			double tr;			// �����X�e�b�v�T�C�Y
			double *dpar;
			int ierr;
			int kd[2];
			int ipar[128];
		} set[2];

		map<double, double>wxt;				// �p���x���ԕω�(first: ����[s], second:�p���x[rad/s])
		map<double, double>vxt;				// ���x���ԕω��@(first: ����[s], second:���x[m/s])
		map<double, vector<double>>ft;		// �׏d���ԕω��@(first: ����[s], second:�׏d[N]/[Nm])

		struct Load {
			//bool is_change;		// ���x������؂�ւ��邩�ǂ����i���ω��X�e�b�v2�ȏ�j
			struct Param {
				double t;			// [s]:		���x�������؂�ւ�鎞��
				Vector3d F;		// [N]:		t�̏u�Ԃ̊O���׏d
				Vector3d M;		// [Nm]:		t�̏u�Ԃ̊O���׏d
			};
			vector<Param> param;
			int i0;				// ���݂̋��
		} lod;
	}dyn;

public:
	static void init_stt(const BS_FileIn::Static & stt, int ballnum);

private:
	static void Stt_Eq0(int * m, int * n, double * x, double * f);
	static void Stt_Eq1(int * m, int * n, double * x, double * f);
	static void Stt_Eq2(int * m, int * n, double * x, double * f);

public:
	static int Stt_solve(int i, double*x);

public:
	static void init_dyn(const BS_FileIn::Dynamic & dyn, int ballnum);
	   
private:
	static void Dyn_Eq0(int * n, double * t, double * y, double * dydt);
	static void Dyn_Eq1(int *n, double *t0, double *y, double *dydt);
	static void Dyn_void(void);
	//static void get_dwdt(double t, double &dwdt);
	static void get_Load(double t, map<double, vector<double>>ft, double *F, double *T);

	//static double get_dvdt(double t, map<double, double> vxt);
	static double get_dvdt(double t, const map<double, double>& vxt);

public:
	static void Dyn_solve(double * x0, double * t, int i);


};
