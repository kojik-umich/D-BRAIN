#pragma once
#define _MAX_CONTACT_ 4
#include "Rigid.h"
#include "B4P_In.h"
#include <Eigen\Dense>
using Eigen::Vector3d;
using Eigen::ArrayXd;

// 抽象保持器クラス．
class bal_Cage : public Rigid {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
public:
	void Nothing(void);		// Rigidの具象クラスであることを示すため，意味のない関数を定義する．
	virtual void init(const B4P_In&FI)=0;
	virtual int get_ContactPoint(const Vector3d&BL_x, double BL_r, int np, Vector3d*x, double*k, int *ptt)=0;	// 玉用保持器として，必ず「玉との接触位置」メソッドを実装すること．これがない具象保持器は認めないものとする．
};


// 冠型保持器クラス．
class bal_SnapCage : public bal_Cage {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
	enum ContactPattern { exception, face, edgeo, edgei, aperture, cornero, corneri };

	ContactPattern how_Contact(const Vector3d & BL_x_, int np, Vector3d&dir, double&zo_th, double&zi_th, double&cos_th, double&sin_th, double&cos_gamma, double&sin_gamma);
	int where_Contact(ContactPattern c, int np, const Vector3d&dir, double zo_th, double zi_th, double cos_th, double sin_th, double cos_gamma, double sin_gamma, Vector3d*x_, double*k);

public:
	void init(const B4P_In&FI);
	int get_ContactPoint(const Vector3d & BL_x, double BL_r, int np, Vector3d*x, double*k, int *ptt);

public:
	struct Pocket {
		Vector3d x;			// ポケット球の中心座標 [m]．（保持器幾何中心座標系）
		double x_norm;		// 上記のノルム [m]．
		Vector3d er;		// 上記の方向ベクトル(=Z)
		Vector3d eth;		// 周方向ベクトル(=Y)
		double R;			// ポケットR [m]．
		double h0;			// ポケットR中心から開口部までの高さ [m]．
		double kface;		// 面接触剛性[N/m]．
		double kci;			// 角接触剛性[N/m].
		double kco;			// 角接触剛性[N/m].
		double kopen;		// 開口部剛性[N/m].
		double kedgei;		// エッジ剛性[N/m].←ボールと保持器の接触剛性なので本来ならばBallCagePairが所持するべき数値？？
		double kedgeo;		// エッジ剛性[N/m].←ボールと保持器の接触剛性なので本来ならばBallCagePairが所持するべき数値？？

		double ropen;
		Vector3d corni0;	// 内径側角座標 [m]（11時寄り）（保持器幾何中心座標系）
		Vector3d corni1;	// 内径側角座標 [m]（ 1時寄り）（保持器幾何中心座標系）
		Vector3d corno0;	// 外径側角座標 [m]（11時寄り）（保持器幾何中心座標系）
		Vector3d corno1;	// 外径側角座標 [m]（ 1時寄り）（保持器幾何中心座標系）
		double cos_alp;		// 開口部開き角 α の cosine.
		double cos_gammai;	// 開口部内径側上限 γ の cosine．負の値になるはず．
		double cos_gammao;	// 開口部外径側上限 γ の cosine．正の値になるはず．
		double zimin;		// 内径側エッジの中心からの最小垂直距離 [m]．
		double zimax;		// 内径側エッジの中心からの最大垂直距離 [m]．
		double ziave;		// 内径側エッジの中心からの平均垂直距離 [m]．
		double ziamp;		// 内径側エッジの中心からの垂直距離振幅 [m]．
		double zomin;		// 外径側エッジの中心からの最小垂直距離 [m]．
		double zomax;		// 外径側エッジの中心からの最大垂直距離 [m]．
		double zoave;		// 外径側エッジの中心からの平均垂直距離 [m]．
		double zoamp;		// 外径側エッジの中心からの垂直距離振幅 [m]．
	} *PK;				// 各ポケットの状態を表す構造体．

	int      Z;			// ポケット数．
	double   pcd;		// Pitch Center Diameter [m]
	double   ro;		// 外径の半径 [m]
	double   ri;		// 内径の半径 [m]
	Vector3d xc;		// 保持器重心(x)から見た幾何中心(c)位置ベクトル [m]（保持器座標系）※幾何中心：ポケット中心の成す平面の中心．

public:
	bal_SnapCage();
	~bal_SnapCage();
};

