#pragma once
#include <Eigen\Dense>
#include "Ball.h"
#include "Tribology.h"
#include "BS_In.h"
#include "BS_Out.h"
#include "BS_Cylinder.h"
#include "BS_Nut.h"
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Matrix3d;


class BS_BallCylinderPair {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:
	Ball       *BL;		// 登録されている玉．
	BS_Cylinder*CY;		// 登録されている螺旋円筒部材．
	int         iSP;		// 条番号．
	double       E;		// 等価ヤング率 [Pa]
	double       m;		// 換算質量 [kg]

	// 接触計算に使うトライボロジー式
	struct Tribologies {
		Tribology::Hertz			*HZ;	// ヘルツ接触式
		Tribology::RollingResistance*RR;	// 転がり粘性抵抗式
		Tribology::FilmThickness	*FT;	// 油膜厚さ式
		Tribology::Traction			*TR;	// 油膜トラクション力の式
		Tribology::Coulomb			*CL;	// クーロン摩擦式
		double						cs;		// クーロン摩擦の立ち上がり傾き．
		Tribology::Hysteresis		*HY;	// ヒステリシス損失
		double						fh;		// ヒステリシス係数
	}TB;

	// 潤滑状態
	struct Lubrication {
		double eta;		// 粘度 [Pa*s]
		double beta;	// 温度粘度係数 [K^-1]
		double k;		// 油熱伝導率 [W/(K*m)]
		double alpha;	// 圧力粘度係数 [Pa^-1]
		double lm;		// メニスカス長さ [m]
	}LB;

	// 各溝の物性値
	struct Groove {
		double	sigma;		// 表面粗さ [m]
		int		n;			// スライス数
		double  *r;			// 接触圧力の分布比 [-]
		Vector2d*xy;		// スライス座標比 [-] xyは単位円の内側の値．0成分:a方向．1成分:b方向．
		Vector3d*ps;		// スライス座標 [m]（本来はローカル変数として扱いたいが，毎回newするとコストがかかるためメンバとして使用）
		double	mu;			// 摩擦係数 [-]
		double	zeta;		// 減衰比 [-]
	} GV[2];
	
	// 出力保存用．
	struct Save {
		Vector3d eta;			// 螺旋位相角[rad], eta座標[m], zeta座標[m] の配列
		// 各溝の一時変数
		struct Groove {
			Vector3d Fn;		// 垂直荷重 [N]（慣性座標系）
			Vector3d Fs;		// 滑り摩擦（スライスの総和） [N]（慣性座標系）
			Vector3d us;		// 滑り速度 [m/s]（慣性座標系）
			Vector3d ur;		// 転がり速度 [m/s]（慣性座標系）
			double dx;			// 弾性接近量 [m]
			double phi;			// 接触角 [rad]
			double a;			// 接触楕円長径 [m]
			double b;			// 接触楕円短径 [m]
			double h;			// 油膜厚さ [m]
			double fratio;		// 油膜荷重支持 [-]
			double lambda;		// ラムダ値 [-]
			Vector3d p;			// 接触楕円中央点 [m]（慣性座標系）
			double F;			// 全ての力の総和 [N]
			Vector2d Fbc;		// 垂直荷重 [N]（溝直行断面座標系）
			double Pmax;		// 最大面圧[Pa]
			struct Slice {
				double f_arr;		// ヘルツ接触荷重[N](スカラー)
				Vector3d ps;		// スライス片中心(慣性座標系)
				Vector3d us;		// 滑り速度[m/s](慣性座標系)
				Vector3d fs;		// 滑り摩擦力[N](慣性座標系)
				Vector3d ts;		// 滑り摩擦トルク[Nm](慣性座標系)
				double mu_cl;		// クーロン摩擦係数[-]
				double mu_tr;		// トラクション係数[-]
			} *SL;					// 各スライス片の結果（ポインタ配列，長さ msmax）
		} GV[2];
	} SV;

public:
	void link(Ball * BL, BS_Cylinder * CY, int is);

	bool how_Contact(int i, const Vector2d & bl_x, Vector2d & e, double & dx);
	double calc_Hertz(int ig, const Vector2d & bl_eta, const Vector2d & e, double dx, double & Rx, double & Ry, Vector2d & p, double & cos_alp, double & sin_alp, double & a, double & b, double & k, Vector2d & rho);
	double calc_DynamicHertz(int ig, const Vector2d & bl_eta, const Vector2d & bl_etav, const Vector2d & e, double dx, double & Rx, double & Ry, Vector2d & p, double & cos_alp, double & sin_alp, double & a, double & b, Vector2d & rho);
	void get_F0(Vector2d & Fbc, Vector3d & Fcb, Vector3d & Tcb);
	void save_F0(int i, const Vector3d & p, double F, const Vector3d & eta, double a, const Vector2d & Fbc);
	void get_F1(bool v, Vector2d & Fbc, Vector3d & Fcb, Vector3d & Tcb);
	void get_F2(Vector3d & vF, Vector3d & vT);
	void get_FT(Vector3d & Fbc, Vector3d & Tbc, Vector3d & Tcb, Vector3d & Fs, Vector3d & Ts);
	Vector3d get_ur(const Vector3d & p);
	Vector3d get_us(const Vector3d & p);
	Vector3d get_eta(void);
	void calc_Sliding(int ig, double th, const Vector3d & p, const Vector3d & xai, double a, double b, double Pmax, double ur_norm, double F_norm, double fratio, const Vector2d & rho, Vector3d & Fs, Vector3d & Ts, Vector3d & Tcs);
	void init(const BS_In::BallCylinderPair & BCP, const BS_In::Tribology & tribology, const BS_In::Oil & oil);
	void init_Tribology(const BS_In::Tribology & tribology);
	void init_Oil(const BS_In::Oil & oil);
	void save(int ig, const Vector3d & eta, const Vector3d & Fn, const Vector3d & Fs, const Vector3d & us, const Vector3d & ur, double dx, double cos_alp, double sin_alp, double a, double b, double h, const Vector3d & p, double lambda, double fratio, double Pmax);
	Vector3d get_e(void);
	void save(BS_Out::BallCylinderPair & OUT);
	void save_Slice(int i, int j, double farr, const Vector3d& Fs_, const Vector3d& Ts_, double mu_cl, double mu_tr, const Vector3d& us, const Vector3d& ps);
	void init_Sliceparam(int i);
	void write_slice(int ig, Matrix3d xyz2XYZ, BS_Out::BallCylinderPair&BRP);
	BS_BallCylinderPair();
	~BS_BallCylinderPair();
};


class BS_BallNutPair : public BS_BallCylinderPair {

public:
	double      th0;	// （主に静解析で利用する）ナットから見た玉の位相角．
	void set_eta0(const Vector2d&eta);
	Vector2d get_etavector0(const Vector3d & x);
};

class BS_BallShaftPair : public BS_BallCylinderPair {
public:
	double      thp;	// （主に静解析で利用する）ナットから見た玉の位相角．

	void set_etap(const Vector2d&eta);
	Vector2d get_etavectorp(const Vector3d & x);
};


//public:
//	enum Contact {
//		NC,		// 接していない状態．NotContact．
//		IVR,	// 非弾性流体潤滑状態．現在実装していない．
//		EHL		// 弾性流体潤滑状態．
//	};
