#pragma once
#include <iostream>
#include <Eigen\Dense>
#include "Ball.h"
#include "Tribology.h"
#include "bal_Cage.h"
#include "B4P_In.h"
#include "B4P_Out.h"
using Eigen::Vector3d;
using Eigen::ColPivHouseholderQR;

// ボール-保持器間のループ不変量を保存しておくクラス．CG=cage．
class B4P_BallCagePair{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
private:
	Tribology::DampingForce*DF;

	Ball       *BL;
	bal_Cage   *CG;
	int         np;		// ボールが対応する保持器ポケット番号．
	double      mu;		// 保持器-ボールのクーロン摩擦係数．(typical:0.1~0.3)
	double      zeta;	// 保持器-ボールのダンピング係数．(typical:~0.2)
	double      m;		// 換算質量．
	Vector3d    p[_MAX_CONTACT_];	// 接触点位置（慣性座標系）
	Vector3d	Fn[_MAX_CONTACT_];	// 接触点荷重（慣性座標系）
	Vector3d	Fs[_MAX_CONTACT_];	// 接触点摩擦（慣性座標系）
	double		dx[_MAX_CONTACT_];	// 接近量[m]
	int			ptt[_MAX_CONTACT_];	// 接触パターン（本当はenumで管理すべきだが，保持器によって定義が異なるため，intで管理）


public:
	void init(const B4P_In&FI);
	void link(Ball*BL, bal_Cage*CG, int np);
	virtual void calc_force(Vector3d&Fbc,Vector3d&Tbc,Vector3d&Fcb,Vector3d&Tcb);
	void save(B4P_Out::BallCagePair&BCP);
	Vector3d get_us(const Vector3d&p);
};

