#pragma once
#include <Eigen\Dense>
#include <iostream>
#include "Ball.h"
#include "B4P_Ring.h"
#include "Tribology.h"
#include "B4P_In.h"
#include "B4P_Out.h"
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::ColPivHouseholderQR;
using Eigen::AngleAxisd;
using namespace std;

// ボール-リング間のループ不変量を保存しておくクラス．
class B4P_BallRingPair {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	Ball    *BL;
	B4P_Ring*RG;
	double sigma;	// リング-ボールの合成粗さ．
	double mu;		// リング-ボールのクーロン摩擦係数．(typical:~0.1)
	double zeta;	// リング-ボールのダンピング係数．(typical:~0.2)
	double K;		// リング-ボールの非線形ばね定数．
	double E;		// リング-ボールの等価ヤング率．
	double lm;		// メニスカス長さ．
	double m;		// 換算質量．
	double lb_alpha;
	double lb_beta;
	double lb_eta;
	double lb_k;
	double clmb_s;	// クーロン摩擦のタンジェントカーブの傾き
	double fh;		// ヒステリシス損失係数

	Vector3d get_ur(const Vector3d&p);
	Vector3d get_us(const Vector3d&p);



	// 外部書き込み用変数群
	Vector3d thXZ;
	struct Slice {
		// スライス片結果
		double f_arr;		// ヘルツ接触荷重[N](スカラー)
		Vector3d ps;		// スライス片中心(慣性座標系)
		Vector3d us;		// 滑り速度[m/s](慣性座標系)
		Vector3d fs;		// 滑り摩擦力[N](慣性座標系)
		Vector3d ts;		// 滑り摩擦トルク[Nm](慣性座標系)
		double mu_cl;		// クーロン摩擦係数[-]
		double mu_tr;		// トラクション係数[-]
	};


	struct Groove {
		Vector3d p;			// 接触点位置（慣性座標系）
		Vector3d Fn;		// 接触荷重（慣性座標系）
		Vector3d us;
		Vector3d ur;
		Vector3d Tr;		// 転がりトルク（慣性座標系）
		Vector3d Fs;
		Vector3d Ts;
		Vector3d Fb;
		Vector3d Tb;
		Vector3d Ti;
		double a;
		double b;
		double cos_alp;
		double sin_alp;
		double Pmax;		// 最大面圧[Pa]
		double lambda;		// ラムダ値[-]
		double fratio;		// 油膜接触割合[-]

		int msmax;			// スライス分割数[-]
		double dx;			// 玉接近量[m]

		Vector2d p_;		// 接触点位置（溝直行断面座標系）
		Vector3d Fn_;		// 接触荷重（溝直行断面座標系）
		Vector3d us_;
		Vector3d ur_;
		Vector3d Tr_;		// 転がりトルク（溝直行断面座標系）
		Vector3d Fs_;
		Vector3d Ts_;
		Vector3d Fb_;
		Vector3d Tb_;
		Vector3d Ti_;
		Slice *SL;
		// 以下の変数はメモリ確保の関係上，"例外的に"計算途中でも用いる
		Vector3d *ps;		// 各スライスの中心位置[m]
		double *ratio_slice;	// 各スライスの接触荷重分布(スカラー)[-]
	}GV[2];
public:
	void link(Ball*BL, B4P_Ring*RG);
	bool how_Contact(int i, const Vector3d&bl_x,  Vector3d&er, Vector3d&eg, double&dx);
	void calc_Hertz(int i, const Vector3d & bl_x, const Vector3d & er, const Vector3d & eg, double dx, double & Rx, double & Ry, Vector3d & p, double & cos_alp, double & sin_alp, double & a, double & b, double & k);
	double calc_DynamicHertz(int i, const Vector3d & bl_x, const Vector3d & bl_v, const Vector3d & er, const Vector3d & eg, double dx, double & Rx, double & Ry, Vector3d & p, double & cos_alp, double & sin_alp, double & a, double & b);
	void calc_Sliding(int i, const Vector3d & p, double a, double Pmax, double F_norm, double fratio, Vector3d & Fs, Vector3d & Tbs, Vector3d & Tis);
	void save(bool c, int i, double cos_alp, double sin_alp, double dx, double a, double b, const Vector3d & p, const Vector3d & Fn, double Pmax, const Vector3d & us, const Vector3d & ur, const Vector3d & Tr, double lambda, double fratio, const Vector3d & Fs, const Vector3d & Ts, const Vector3d & Fb, const Vector3d & Tb, const Vector3d & Ti);
	void save_Slice(int i, int j, double farr, const Vector3d& Fs_, const Vector3d& Ts_, double mu_cl, double mu_tr, const Vector3d& us, const Vector3d& ps);
	void init_Sliceparam(int i);
	virtual void calc_force(Vector3d&Fbi, Vector3d&Nbi, Vector3d&Fib, Vector3d&Nib);
	void write(B4P_Out::BallRingPair&BXP);
	void write_slice(int ig, Matrix3d xyz2XYZ, B4P_Out::BallRingPair&GV);
	void init_Tribology(const B4P_In::Tribology & FI);
	void init_Lubrication(const B4P_In::Lubrication & FI);
	void init_Slice(int msmax);
	virtual void init(const B4P_In&FI) = 0;
	virtual double ContactAngle(double cos_alp, double sin_alp) = 0;
	void get_us_ur(const Vector3d&p, Vector3d&us, Vector3d&ur);
	void get_us_ur2(const Vector3d&p, const Vector3d&er, Vector3d&us, Vector3d&ur);
	void get_Fstf(Vector3d& Fbi, Vector3d& Fib, Vector3d& Tib);
	B4P_BallRingPair();
	~B4P_BallRingPair();

	// 動的に切り替わるトライボロジークラス．
	Tribology::Hertz *HZ;
	Tribology::RollingResistance *RR;
	Tribology::FilmThickness*FT;
	Tribology::Traction*TR;
	Tribology::Coulomb*CL;
	Tribology::DampingForce*DF;
	Tribology::Hysteresis *HY;
};

class B4P_BallOuterRingPair : public B4P_BallRingPair {
public:
	void init(const B4P_In&FI);
	double ContactAngle(double cos_alp, double sin_alp);
	void set_eta_stf(const Vector3d & eta);
};

class B4P_BallInnerRingPair : public B4P_BallRingPair {
public:
	void init(const B4P_In&FI);
	double ContactAngle(double cos_alp, double sin_alp);

};



//VectorXd make_SliceParam(void);
//double get_c(double k);
		//double*fn_slice;	// 各スライスの接触荷重(スカラー)[N]
		//Vector3d*Ts_slice;	// 各スライスのすべり摩擦モーメント[N]
		//Vector3d *Fs_slice;	// 各スライスの滑り摩擦[N]
		//double *mutr_slice;	// 各スライスのトラクション係数[-]
		//double *mucl_slice;	// 各スライスのトラクション係数[-]
		//Vector3d *us_slice;	// 各スライスの滑り速度[m/s]
