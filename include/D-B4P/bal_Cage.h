#pragma once
#define _MAX_CONTACT_ 4
#include "Rigid.h"
#include "B4P_In.h"
#include <Eigen\Dense>
using Eigen::Vector3d;
using Eigen::ArrayXd;

// ���ەێ���N���X�D
class bal_Cage : public Rigid {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
public:
	void Nothing(void);		// Rigid�̋�ۃN���X�ł��邱�Ƃ��������߁C�Ӗ��̂Ȃ��֐����`����D
	virtual void init(const B4P_In&FI)=0;
	virtual int get_ContactPoint(const Vector3d&BL_x, double BL_r, int np, Vector3d*x, double*k, int *ptt)=0;	// �ʗp�ێ���Ƃ��āC�K���u�ʂƂ̐ڐG�ʒu�v���\�b�h���������邱�ƁD���ꂪ�Ȃ���ەێ���͔F�߂Ȃ����̂Ƃ���D
};


// ���^�ێ���N���X�D
class bal_SnapCage : public bal_Cage {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
	enum ContactPattern { exception, face, edgeo, edgei, aperture, cornero, corneri };

	ContactPattern how_Contact(const Vector3d & BL_x_, int np, Vector3d&dir, double&zo_th, double&zi_th, double&cos_th, double&sin_th, double&cos_gamma, double&sin_gamma);
	int where_Contact(ContactPattern c, int np, const Vector3d&dir, double zo_th, double zi_th, double cos_th, double sin_th, double cos_gamma, double sin_gamma, Vector3d*x_, double*k);

public:
	void init(const B4P_In&FI);
	int get_ContactPoint(const Vector3d & BL_x, double BL_r, int np, Vector3d*x, double*k, int *ptt);

public:
	struct Pocket {
		Vector3d x;			// �|�P�b�g���̒��S���W [m]�D�i�ێ���􉽒��S���W�n�j
		double x_norm;		// ��L�̃m���� [m]�D
		Vector3d er;		// ��L�̕����x�N�g��(=Z)
		Vector3d eth;		// �������x�N�g��(=Y)
		double R;			// �|�P�b�gR [m]�D
		double h0;			// �|�P�b�gR���S����J�����܂ł̍��� [m]�D
		double kface;		// �ʐڐG����[N/m]�D
		double kci;			// �p�ڐG����[N/m].
		double kco;			// �p�ڐG����[N/m].
		double kopen;		// �J��������[N/m].
		double kedgei;		// �G�b�W����[N/m].���{�[���ƕێ���̐ڐG�����Ȃ̂Ŗ{���Ȃ��BallCagePair����������ׂ����l�H�H
		double kedgeo;		// �G�b�W����[N/m].���{�[���ƕێ���̐ڐG�����Ȃ̂Ŗ{���Ȃ��BallCagePair����������ׂ����l�H�H

		double ropen;
		Vector3d corni0;	// ���a���p���W [m]�i11�����j�i�ێ���􉽒��S���W�n�j
		Vector3d corni1;	// ���a���p���W [m]�i 1�����j�i�ێ���􉽒��S���W�n�j
		Vector3d corno0;	// �O�a���p���W [m]�i11�����j�i�ێ���􉽒��S���W�n�j
		Vector3d corno1;	// �O�a���p���W [m]�i 1�����j�i�ێ���􉽒��S���W�n�j
		double cos_alp;		// �J�����J���p �� �� cosine.
		double cos_gammai;	// �J�������a����� �� �� cosine�D���̒l�ɂȂ�͂��D
		double cos_gammao;	// �J�����O�a����� �� �� cosine�D���̒l�ɂȂ�͂��D
		double zimin;		// ���a���G�b�W�̒��S����̍ŏ��������� [m]�D
		double zimax;		// ���a���G�b�W�̒��S����̍ő吂������ [m]�D
		double ziave;		// ���a���G�b�W�̒��S����̕��ϐ������� [m]�D
		double ziamp;		// ���a���G�b�W�̒��S����̐��������U�� [m]�D
		double zomin;		// �O�a���G�b�W�̒��S����̍ŏ��������� [m]�D
		double zomax;		// �O�a���G�b�W�̒��S����̍ő吂������ [m]�D
		double zoave;		// �O�a���G�b�W�̒��S����̕��ϐ������� [m]�D
		double zoamp;		// �O�a���G�b�W�̒��S����̐��������U�� [m]�D
	} *PK;				// �e�|�P�b�g�̏�Ԃ�\���\���́D

	int      Z;			// �|�P�b�g���D
	double   pcd;		// Pitch Center Diameter [m]
	double   ro;		// �O�a�̔��a [m]
	double   ri;		// ���a�̔��a [m]
	Vector3d xc;		// �ێ���d�S(x)���猩���􉽒��S(c)�ʒu�x�N�g�� [m]�i�ێ�����W�n�j���􉽒��S�F�|�P�b�g���S�̐������ʂ̒��S�D

public:
	bal_SnapCage();
	~bal_SnapCage();
};

