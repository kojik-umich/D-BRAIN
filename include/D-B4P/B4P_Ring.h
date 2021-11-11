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

/// <summary>4点接触玉軸受に使われるリングクラス．内外輪どちらも含む．</summary>
/// <remarks>中心の違う，２つのトーラス溝をもつリング．トーラスは特にクラス化などはさせていない．必要あれば今後の修正で独立化させても大丈夫です．</remarks>
class B4P_Ring : public Rigid {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
public:
	double de;			// 外径 [m]
	double di;			// 内径 [m]
	double pcd;			// ピッチ円直径 [m]
	double sigmap;		// 粗さrms(新規追加)
	struct Groove {
		double Rx;			// 内外輪中心と溝R中心のx方向ずれ量 （0:負，1:正）[m]
		double r;			// 溝径 [m]
		double r_inv;		// 溝径の曲率 [1/m]
		double Rr;			// トーラス円半径（内外輪中心から溝中心の距離，溝R中心pcd÷2） [m]
		double h;			// 溝肩高さ[m](新規追加)
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


