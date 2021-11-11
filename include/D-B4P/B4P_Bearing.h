#pragma once
#include <stdio.h>
#include <Eigen\Dense>
#include "B4P_In.h"
#include "B4P_Out.h"
#include "Numeric.h"
#include "Ball.h"
#include "B4P_Ring.h"
#include "B4P_BallRingPair.h"
#include "bal_Cage.h"
#include "B4P_BallCagePair.h"
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix3d;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::ArrayXd;
using Eigen::ColPivHouseholderQR;
using Eigen::Map;

/// <summary>4点接触玉軸受</summary>
/// <remarks>各部材（玉・内外輪・保持器）を取りまとめるクラス．現在の軸受内部の状態を計算できる．</remarks>
class B4P_Bearing {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:
	Ball	 *BL;
	B4P_OuterRing OR;
	B4P_InnerRing IR;
	bal_Cage *CG;
	B4P_BallOuterRingPair *BOP;
	B4P_BallInnerRingPair *BIP;
	B4P_BallCagePair *BCP;
	Vector3d F_load;	// 外部荷重 [N]．
	Vector3d T_load;	// 外部トルク [Nm]． 
	double   pcd;		// Pitch Center Diameter [m]
	double *azimuth0;	// 方位角 [rad]

public:
	int    Z;			// 球数．
	int    nX;			// 変位ベクトル長さ．（動解析用）
	void init(const B4P_In & FI, double * x, double * ax);
	void allocate_Pair(int Z);
	void set_y(const double*y);
	void get_y(double*y);
	void get_dydt(double*dydt);
	void save(B4P_Out&OUT);
	void set_Xstf(double* Xstf);
	void get_Fstf(double* Fstf);
	void get_Xstf(double* Xstf);


private:
	void init_parts(const B4P_In & FI, double*Cage_rmg0);
	void init_position(double balldia, int ballnum, double ballpcd, double cos_alp0, double omegair, double omegaor, const double*rmg0);

public:
	B4P_Bearing(void);
	~B4P_Bearing(void);
};




//double get_fn_ball();
// Vector3d get_myz_ir();
// double get_fn_ir();
// double get_fx_ball();
// double get_fx_ir();
