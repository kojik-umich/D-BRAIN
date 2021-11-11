#pragma once
#include <Eigen\Dense>
#include "Spiral.h"
#include "BS_In.h"

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Matrix3d;

class BS_Spiral {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
private:
	Spiral SP;		// �����N���X

public:
	struct Groove {
		Vector2d eta;	// �a���S�̌a�����Y�� [m]�D����+�D�i�b�g���C�V���t�g���D
		double r;		// �a���a [m]�D
		double r_inv;	// �a�ȗ� [1/m]�D
		double sigma;	// �\�ʑe�� [m]
		double R;		// �ȗ�R [m]�D���[�v�s�ϗʂȂ̂Ń����o�ɕۑ��D
	} GV[2];

public:
	// �Ϗ����邾���ő��x�������o��͉̂������̂ŁC�Ȃ�ׂ��x���Ȃ�Ȃ��悤�ɃC�����C���W�J������D
	inline Vector3d to_eta(const Vector3d&xyz)	{ return this->SP.to_eta2(xyz); };
	inline Vector3d to_xyz(const Vector3d&eta)	{ return this->SP.to_xyz(eta); };
	inline Matrix3d get_xyz2eta(double theta)	{ return this->SP.get_xyz2eta(theta); };
	inline double   get_nd(void)				{ return this->SP.get_nd(); };
	inline double   get_r(void)					{ return this->SP.get_r(); };

	void Nothing(void);
	void init(const BS_In::Cylinder::Spiral & spiral);
	Vector2d get_rho(double cos_alp, int i);
};




