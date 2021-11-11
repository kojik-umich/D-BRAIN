/*******************************************************************************
!								"BS_Nut.cpp"
!													2019/03/22	[Core-T]	楠崎
!
!
!*******************************************************************************/
#include "BS_Nut.h"
using Eigen::VectorXd;


BS_Nut::BS_Nut() {
	this->CY = NULL;
	return;
}

BS_Nut::~BS_Nut() {
	if (this->CY != NULL)
		delete[] this->CY;
	return;
}

// 以下，シングルナット．
void BS_SingleNut::Nothing(void) {
	return;
}

void BS_Nut::allocate(const std::vector<BS_In::Cylinder>&cylinders) {

	this->nCY = cylinders.size();
	this->CY = new BS_Cylinder[this->nCY];

	for (int i = 0; i < this->nCY; i++)
		this->CY[i].allocate(cylinders[i].spiral.size());
	return;
}

void BS_SingleNut::init(const std::vector<BS_In::Cylinder> & nut, double w0) {
	bool v_const[3], w_const[3];
	for (int i = 0; i < 3; i++) {
		v_const[i] = true;
		w_const[i] = true;
	}

	this->set_const(v_const, w_const);


	Vector3d x, v, w; Quaterniond q;
	x = Vector3d::Zero();
	v = Vector3d::Zero();
	w = Vector3d::Zero();
	q = Quaterniond::Identity();

	this->set_y_(x, v, q, w);

	this->CY[0].init(nut[0], v_const, w_const);
	return;
}

void BS_SingleNut::set_y_(const Vector3d & x, const Vector3d & v, const Quaterniond & q, const Vector3d & w) {
	this->set_y(x, v, q, w);
	this->CY[0].set_y(x, v, q, w);
	return;
}

// 以下，ダブルナット．
void BS_DoubleNut::Nothing(void) {
	return;
}



void BS_DoubleNut::init(const std::vector<BS_In::Cylinder> & nut, double w0) {
	bool v_const[3], w_const[3];
	for (int i = 0; i < 3; i++) {
		v_const[i] = true;
		w_const[i] = true;
	}
	this->set_const(v_const, w_const);

	for (int i = 0; i < 2; i++)
		this->dx[i] = Vector3d(nut[i].x0, 0.0, 0.0);

	Vector3d x, v, w; Quaterniond q;
	x = Vector3d::Zero();
	v = Vector3d::Zero();
	w = Vector3d(w0, 0.0, 0.0);
	q = Quaterniond::Identity();

	this->set_y(x, v, q, w);
	this->set_dx();

	for (int i = 0; i < 2; i++)
		this->CY[i].init(nut[i], v_const, w_const);

	return;
}

// ナット中心から左右のナットそれぞれの変位を与えるメソッド．
void BS_DoubleNut::set_dx(void) {

	for (int i = 0; i < 2; i++) {
		Vector3d dx = this->to_inecoord(this->dx[i]);
		this->CY[i].set_param(this->x + dx, this->v, this->q, this->w);
	}
	return;
}

// 値の設定メソッド．無次元量で渡された値を有次元にしてパラメタとして代入する．
void BS_DoubleNut::set_y_(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w) {

	this->set_y(x, v, q, w);
	this->set_dx();

	return;
}






























//#include "BS_Screw.h"
//
//
//void BS_Screw::init(int ns) {
//
//	this->SP = new BS_Spiral[ns];
//}
//
//// ナットの溝直交座標系に変換するメソッド．
//Vector3d BS_Screw::to_etacoord
//(						// out:	[rad],[m]:	0成分：ナット位相角．1成分：eta座標．2成分：zeta座標．
//	const Vector3d&x,	// in : [m]:		変換したい慣性座標系．
//	int i				// in : [-]:		変換したい条番号．
//) {
//	Vector3d xyz = this->to_mycoord(x);
//	Vector3d eta = this->SP[i].to_eta2(xyz);
//	return eta;
//}
//
//// ナットの溝直交座標系に変換するメソッド．
//Vector3d BS_Screw::to_inertialcoord
//(						// out:	[m]:		慣性座標系の座標．
//	const Vector3d&eta,	// in :	[rad],[m]:	変換したい螺旋座標系．
//	int i				// in : [-]:		変換したい条番号．
//) {
//	Vector3d xyz = this->SP[i].to_xyz(eta);
//	Vector3d x = this->to_inecoord(xyz);
//	return x;
//}
// 
//void BS_Screw::set_TransMat
//(
//	double theta,	// in : [rad]:		ナット位相角．
//	int i			// in : [-]:		変換したい条番号．
//)
//{
//	this->SP[i].set_TransMat(theta);
//}
//// 慣性座標系速度から螺旋座標系速度へ変換するメソッド．
//// ※使う前に，必ずset_TransMatで更新してください！
//Vector3d BS_Screw::to_etavelocity
//(					// out:	[m/s]:	螺旋座標系における速度3成分．
//	Vector3d v,		// in : [m/s]:	慣性座標系の変換したい速度．
//	int i			// in : [-]:	変換したい条番号．
//)
//{
//	Vector3d v_ = this->to_myvelocity(v);			// ナット座標系への変換．
//	Vector3d v__ = this->SP[i].to_etavector(v_);	// 螺旋座標系への変換．
//	return v__;
//}
//
//// 慣性座標系速度から螺旋座標系速度へ変換するメソッド．
//// ※使う前に，必ずset_TransMatで更新してください！
//Vector3d BS_Screw::to_inertialvelocity
//(					// out:	[m/s]:	慣性座標系における速度3成分．
//	Vector3d v,		// in : [m/s]:	螺旋座標系の変換したい速度．
//	int i			// in : [-]:	変換したい条番号．
//)
//{
//	Vector3d v_ = this->SP[i].to_xyzvector(v);		// ナット座標系への変換．
//	Vector3d v__ = this->to_inevelocity(v_);		// 慣性座標系への変換．
//	return v__;
//}

