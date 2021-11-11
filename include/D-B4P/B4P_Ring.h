#pragma once
#include "Numeric.h"
#include "Rigid.h"
#include "B4P_In.h"
#include <Eigen\Dense>
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::VectorXd;
using Eigen::Quaterniond;

/// <summary>4�_�ڐG�ʎ���Ɏg���郊���O�N���X�D���O�ւǂ�����܂ށD</summary>
/// <remarks>���S�̈Ⴄ�C�Q�̃g�[���X�a���������O�D�g�[���X�͓��ɃN���X���Ȃǂ͂����Ă��Ȃ��D�K�v����΍���̏C���œƗ��������Ă����v�ł��D</remarks>
class B4P_Ring : public Rigid {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
public:
	double de;			// �O�a [m]
	double di;			// ���a [m]
	double pcd;			// �s�b�`�~���a [m]
	double sigmap;		// �e��rms(�V�K�ǉ�)
	struct Groove {
		double Rx;			// ���O�֒��S�ƍaR���S��x��������� �i0:���C1:���j[m]
		double r;			// �a�a [m]
		double r_inv;		// �a�a�̋ȗ� [1/m]
		double Rr;			// �g�[���X�~���a�i���O�֒��S����a���S�̋����C�aR���Spcd��2�j [m]
		double h;			// �a������[m](�V�K�ǉ�)
	} GV[2];
	   
public:
	void Nothing(void);	
	virtual void init(const B4P_In::Ring&IN, double pcd);
	void calc_slice(const Vector3d&p, double a, int msmax, int i, double bl_x, Vector3d*ps);
	void get_dydt_(const Vector3d & F, const Vector3d & T, double * dydt);
	Vector3d XZ_to_inecoord(const Vector3d & eta);
	Vector3d ine_to_XZvector(const Vector3d & a, const Matrix3d & xyz2XYZ);
	Matrix3d get_xyz2XYZ(double th);
	Vector3d ine_to_XZcoord(const Vector3d & x);
	Vector2d to_cross_sction(const Vector3d & x, double th);
	Matrix3d get_Rth();
	void set_y_(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w);
	Vector4d get_dqdt_(bool wy_const, bool wz_const);
	Vector3d get_dwdt_(const Vector3d&T, bool wy_const, bool wz_const);
	//void set_constparam(const Vector3d&x, const Vector3d&ax);
};

class B4P_OuterRing : public B4P_Ring {
public:
};

class B4P_InnerRing : public B4P_Ring {
public:
	void set_Xstf(double *x);
	void get_Xstf(double* Xstf);
};


