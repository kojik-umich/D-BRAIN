#pragma once
#include <Eigen\Dense>	// �����p�����^��Vector3d�p�D
#include <math.h>		// sin, cos �̌v�Z�p�D
#include "Numeric.h"	// pi �̎g�p�D
using Eigen::Vector3d;
using Eigen::Matrix3d;

class Spiral {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
	double alp;		// �����n�_�p[rad]�i���萔�j
	double l;		// ���[�h�� [m]
	double r;		// �����s�b�`���a [m]
	double l_nd;	// ���������[�h�� [-]
	double r_nd;	// �����������s�b�`���a [-]
	double nd;		// ���������̂��߂̗�(Non Dimension) [m]�D����1�s�b�`���̒����ɓ������D
	double l_2pi;	// ���[�h�� / 2�� �̒l [m]

public:
	void init(double alp,double l,double r);
	Vector3d to_xyz(const Vector3d&eta);
	Vector3d to_eta2(const Vector3d&xyz);
	Matrix3d get_xyz2eta(double theta);
	double   get_nd(void);
	double   get_r(void);

private:
	Vector3d to_eta(const Vector3d&xyz);
};


