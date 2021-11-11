#include "Ball.h"

// ダミー関数．
void Ball::Nothing(void) {
}


// Ball([玉]0玉番号，1密度，2ヤング率，3por，4粗さrms，5玉径)のパラメータを代入
void Ball::init(
	double D,		// 玉径
	double E,		// ヤング率
	double por,		// ポアソン比
	double den,		// 密度
	double rms,		// 粗さrms
	const bool(&v_const)[3],	// 拘束条件（速度）
	const bool(&w_const)[3]		// 拘束条件（角速度）
) {
	this->D				= D;	
	this->E				= E;	
	this->nu			= por;	
	this->rho			= den;	
	this->sigmap		= rms;	
	this->m				= Numeric::pi / 6 * D * D * D * rho;
	this->m_inv			= 1.0 / this->m;
	this->I				= this->m*D*D/10 * Vector3d::Ones();
	this->I_inv			= 1.0 / this->I[0] * Vector3d::Ones();
	this->r				= 0.5 * D;
	this->r_inv			= 2.0 / D;
	this->set_const(v_const, w_const);
	return;
}

// 玉に接近する物体の座標・速度を元に，接触の有無，接近量，接近速度，接近方向を計算
bool Ball::calc_Contact(	// out:	[-]:	接触の有無（接触あり：true，接触なし：false）
	const Vector3d&x,		// in:	[m]:	接近する物体の座標．（慣性座標系）
	const Vector3d&v,		// in:	[m/s]:	接近する物体の表面速度．（慣性座標系）
	double &dx,				// out:	[m]:	接近量（スカラー，干渉していれば正）
	double &dv,				// out:	[m/s]:	接近速度（スカラー，接近量が大きくなる向きが正）
	Vector3d&edir			// out:	[-]:	接近方向単位ベクトル．
) {

	Vector3d Zero = Vector3d::Zero();

	// ボール中心から侵入者への距離と，その方向の単位ベクトル．
	double Distance = (x - this->x).norm();
	if (Distance == 0)
		return false;

	// ボール中心から侵入者に向かうベクトル．ボールにかかる荷重と逆向き．
	edir = (x - this->x) / Distance;

	dx = this->r - Distance;

	// 接触していない場合荷重0で返す．
	if (dx < 0)
		return false;

	// 接近する物体の速度 dv を算出する．（スカラー，接近する向きが正）
	dv = (v - this->v).dot(-edir);

	return true;
}



// 入力されたベクトルの，法線方向成分を除去するメソッド．
Vector3d Ball::remove_Normal(	// out:	[]:		法線成分を除去したベクトル．（慣性座標系）
	const Vector3d&p,			// in:	[m]:	侵入者の座標．（慣性座標系）
	const Vector3d&a			// in:	[]:		加工したいベクトル．（慣性座標系）
) {
	Vector3d e = (p - this->x).normalized();
	double ea = e.dot(a);
	return a - ea * e;
}

void BallBallPair::link(Ball * BL0, Ball * BL1) {
	this->BL[0] = BL0;
	this->BL[1] = BL1;
	return;
}
