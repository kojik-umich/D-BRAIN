/*******************************************************************************
!								"Spiral.cpp"
!													2019/12/25	[Core-T]	楠崎
!	螺旋に関する計算のオブジェクト．
!	直交座標系⇔螺旋座標系の変換が主な役割．
!	軸はx軸を向いているものとする．またx=0での平面での始点角をαとし，右ねじで螺旋は進むものとする．
!*******************************************************************************/
#include "Spiral.h"


void Spiral::init(
	double alp,		// in:	[rad]:	螺旋始点位相角．
	double l, 		// in:	[m]:	リード長さ．
	double r		// in:	[m]:	螺旋半径．
) {
	this->alp = alp;
	this->l   = l;
	this->r   = r;

	this->l_2pi = l / (2 * Numeric::pi);
	this->nd    = sqrt(this->r*this->r + this->l_2pi*this->l_2pi);
	this->l_nd  = this->l_2pi / this->nd;
	this->r_nd  = this->r     / this->nd;

	return;
}

// 部材座標系から螺旋座標系に変換するメソッド．
Vector3d Spiral::to_eta
(						// out:	[rad],[m]	螺旋座標系
	const Vector3d&xyz	// in:	[m]:		部材座標系（軸がx方向を向いています）
) {
	double x      = xyz[0];
	double y      = xyz[1];
	double z      = xyz[2];
	double phi    = x / this->nd / this->l_nd;
	double beta   = phi + this->alp;
	double s_beta = sin(beta);
	double c_beta = cos(beta);
	double y0     = c_beta * y + s_beta * z;
	double z0     =-s_beta * y + c_beta * z;
	double th     = z0 / (y0 + this->l_nd * this->l_2pi / this->r_nd);
	double ze     = -th * this->l_2pi / this->r_nd;
	double et     = this->r_nd * this->nd - y0 - y0 * th*th / 2 + this->l_nd * th * ze;
	Vector3d eta  = Vector3d(th+phi, et, ze);
	return eta;
}

// 螺旋座標系から部材座標系に変換するメソッド．
Vector3d Spiral::to_xyz
(							// out:	[m]:		直交座標系(x,y,z)(物体座標系)
	const Vector3d&eta		// in:	[rad],[m]	螺旋座標系(θ,η,ζ)(物体座標系)
) {
	double th     = eta[0];
	double beta   = th + this->alp;
	double s_beta = sin(beta);
	double c_beta = cos(beta);
	Vector3d q    = Vector3d(this->l_2pi * th, this->r * c_beta, this->r * s_beta);
	Vector3d n    = Vector3d(0.0, -c_beta, -s_beta);
	Vector3d b    = Vector3d(this->r_nd, this->l_nd * s_beta, -this->l_nd * c_beta);
	Vector3d xyz  = q + eta[1] * n + eta[2] * b;
	return xyz;
}

// 部材座標系から螺旋座標系に変換するメソッド．1回繰り返し計算を行うことで，そこそこの精度で変換できる．
Vector3d Spiral::to_eta2
(							// out:	[rad],[m]:	螺旋座標系(θ,η,ζ)(物体座標系)
	const Vector3d&xyz		// in:	[m]:		直交座標系(x,y,z)(物体座標系)
) {
	Vector3d eta0 = this->to_eta(xyz);
	Vector3d xyz0 = this->to_xyz(eta0);
	Vector3d xyz1 = 2 * xyz - xyz0;
	Vector3d eta1 = this->to_eta(xyz1);
	return eta1;
}

// x-y-z座標系のベクトルをxi-eta-zeta座標系のベクトルに変換するための行列取得．
Matrix3d Spiral::get_xyz2eta
(							// out: [-]:	3*3行列
	double theta			// in:	[rad]:	基準となる螺旋座標系の原点の位相角
) {
	double beta   = theta + this->alp;
	double s_beta = sin(beta);
	double c_beta = cos(beta);

	Matrix3d xyz2eta;
	xyz2eta <<
		this->l_nd, -this->r_nd * s_beta, this->r_nd * c_beta,
		0.0, -c_beta, -s_beta,
		this->r_nd, this->l_nd * s_beta, -this->l_nd * c_beta;
	return xyz2eta;
}

double Spiral::get_nd(void) {
	return this->nd;
}

double Spiral::get_r(void) {
	return this->r;
}
