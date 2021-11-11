#pragma once
#include <Eigen\Dense>	// 螺旋パラメタのVector3d用．
#include <math.h>		// sin, cos の計算用．
#include "Numeric.h"	// pi の使用．
using Eigen::Vector3d;
using Eigen::Matrix3d;

class Spiral {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
	double alp;		// 螺旋始点角[rad]（＝定数）
	double l;		// リード長 [m]
	double r;		// 螺旋ピッチ半径 [m]
	double l_nd;	// 無次元リード長 [-]
	double r_nd;	// 無次元螺旋ピッチ半径 [-]
	double nd;		// 無次元化のための量(Non Dimension) [m]．螺旋1ピッチ分の長さに等しい．
	double l_2pi;	// リード長 / 2π の値 [m]

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


