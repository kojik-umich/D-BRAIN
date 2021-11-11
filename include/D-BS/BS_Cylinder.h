#pragma once
#include "Numeric.h"
#include "Rigid.h"
#include "BS_Spiral.h"
#include "BS_In.h"
#include "BS_Out.h"
#include <Eigen\Dense>

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::Quaterniond;

class BS_Cylinder : public Rigid {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:
	int    nSP;		// 条数 [-]
	BS_Spiral*SP;	// 条数が分からないため，newで宣言する．

public:
	double ri;		// 内半径 [m]
	double ro;		// 外半径 [m]
	Vector3d Fs;	// 滑りによってかかる荷重 [N]．主に出力記録用．
	Vector3d Ts;	// 滑りによってかかるトルク [Nm]．主に出力記録用．

public:
	void Nothing() {};
	Vector3d to_etacoord(int i, const Vector3d & x);
	Vector3d to_inertialcoord(int i, const Vector3d & eta);
	Matrix3d get_xyz2eta(int i, double theta);
	Vector3d to_etavelocity(const Vector3d & v, const Matrix3d & xyz2eta);
	Vector3d to_inertialvelocity(const Vector3d & v, const Matrix3d & xyz2eta);
	Vector3d to_etavector(const Vector3d & a, const Matrix3d & xyz2eta);
	Vector3d to_inertialvector(const Vector3d & a, const Matrix3d & xyz2eta);
	void calc_slice(int is, int ig, double th, const Vector3d & p, const Vector3d & xai, double a, double b, const Vector2d * xy, int msmax, const Vector2d & Rm, Vector3d * ps);
	void allocate(int ns);
	void init(const BS_In::Cylinder & cylinder, const bool(&v_const)[3], const bool(&w_const)[3]);
	double get_nd(int is);
	double get_r(int is);
	void linspace(int is, double th0, double th1, int nb, Vector3d * xs);
	void save(BS_Out::Cylinder&OUT);

	BS_Cylinder();
	~BS_Cylinder();
};
