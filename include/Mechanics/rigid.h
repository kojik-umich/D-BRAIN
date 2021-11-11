#pragma once
#include "Numeric.h"
#include <Eigen\Dense>
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Quaterniond;

/// <summary>剛体</summary>
/// <remarks>変位や姿勢，質量を持ち，ベクトルの変換や運動方程式による時間微分値を計算する．</remarks>
class Rigid {
public:
	Vector3d x;		// 変位 [m]（0:x,1:y,2:z）
	Vector3d v;		// 速度 [m/s]
	Quaterniond q;	// 姿勢を表すクォータニオン [-]
	Vector3d w;		// 角速度 [rad/s]

	double   m;		// 質量[kg]
	double   m_inv;	// 質量の逆数．
	Vector3d I;		// 慣性モーメント[kg*m^2]
	Vector3d I_inv;	// 慣性モーメントの逆数
	double   rho;	// 密度 [kg/m^3]
	double   E;		// ヤング率 [Pa]
	double   nu;	// ポアソン比 [-]

	Vector3d ax0;	// 初期の軸方向ベクトル．(=[1,0,0])
	Vector3d ay0;	// 方位角180°の径方向ベクトル．造語．(=[0,1,0])
	Vector3d az0;	// 方位角270°の径方向ベクトル．造語．(=[0,0,1])

	Vector3d F;		// 受けている外部荷重[N]．慣性座標系．
	Vector3d T;		// 受けている外部トルク[Nm]．慣性座標系．

	//bool v_const;	// 速度固定
	//bool w_const;	// 角速度固定


	Eigen::Array<bool, 3, 1>x_const;		//	拘束条件（x, y, z成分）
	Eigen::Array<bool, 3, 1>Rx_const;		//	拘束条件（Rx, Ry, Rz成分）

	static double l;	// 単位長さ[m]
	static double t;	// 単位時間[s]
	static Vector3d g;	// 重力加速度[m/ss]

public:
	Rigid(void);
	~Rigid(void);
	virtual void Nothing(void)=0;	// 抽象クラスであることを示すため，意味のない関数を定義する．
	void set_param(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w);
	void get_param(Vector3d&x, Vector3d&v, Quaterniond&q, Vector3d&w);
	void set_y(const Vector3d & x, const Vector3d & v, const Quaterniond & q, const Vector3d & w);
	void get_y(Vector3d & x, Vector3d & v, Quaterniond & q, Vector3d & w);
	Vector3d get_dxdt(void);
	Vector3d get_dvdt(const Vector3d&F);
	Quaterniond get_dqdt(void);
	Vector3d get_dwdt(const Vector3d&T);
	void get_dydt(const Vector3d&F, const Vector3d&T, double*dydt);
	void set_const(const bool (&v_const)[3] ,const bool (&w_const)[3]);
	Vector3d to_mycoord(const Vector3d&x);
	Vector3d to_inecoord(const Vector3d&x);
	Vector3d to_myvelocity(const Vector3d&v);
	Vector3d to_inevelocity(const Vector3d&V);
	Vector3d to_myvector(const Vector3d&x);
	Vector3d to_inevector(const Vector3d&x);
	Vector3d surface_velocity(const Vector3d&x);
	Vector3d calc_Torque(const Vector3d&x, const Vector3d&F);
	Vector3d calc_TorqueDirection(const Vector3d & x, const Vector3d & u);
	Vector3d get_ax(void);
	Vector3d get_ay(void);
	Vector3d get_az(void);
	Vector3d get_mg(void);
	Vector3d remove_ax(const Vector3d&x);
	Vector3d remove_vector(const Vector3d&x);
	void set_ax(const Vector3d&ax);
	void set_FT(const Vector3d&F, const Vector3d&T);
	void set_mI(double m, const Vector3d&I);
	void save(double * x, double * v, double * q, double * w, double * ax, double * F, double * T);
};

