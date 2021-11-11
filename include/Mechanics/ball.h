#pragma once
#include <Eigen\Dense>
#include<iostream> 
#include<algorithm> 
#include "Numeric.h"
#include "Rigid.h"
using Eigen::Vector3d;
using Eigen::Quaterniond;
using std::max; 

class Ball : public Rigid {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:
	double   D;				// ���a[m]
	double   r;				// ���a[m]
	double   r_inv;			// �ȗ�[1/m]
	double	 sigmap;		// �e��rms

public:
	void Nothing(void);	// ��ۃN���X�ł��邱�Ƃ��������߁C�Ӗ��̂Ȃ��֐����`����D
	void init(double D, double E, double por, double den, double rms, const bool(&v_const)[3], const bool(&w_const)[3]);
	bool calc_Contact(const Vector3d&x, const Vector3d&v, double &dx, double &dv, Vector3d &edir);
	Vector3d remove_Normal(const Vector3d&x, const Vector3d&a);

};


class BallBallPair {
public:
	Ball *BL[2];
	void link(Ball*BL0, Ball*BL1);

};

