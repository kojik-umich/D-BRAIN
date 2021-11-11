/*******************************************************************************
!								"BS_Screw.cpp"
!													2020/01/08	[Core-T]	楠崎
!	ねじ部材を表すオブジェクト．
!	これを継承してナット・シャフトを作成する．
!	溝の条数分だけ new で作成する．
!
!*******************************************************************************/
#include "BS_Cylinder.h"




// 座標を慣性座標系からナットの溝直交座標系へ変換
Vector3d BS_Cylinder::to_etacoord
(						// out:	[rad],[m]:	0成分：ナット位相角．1成分：eta座標．2成分：zeta座標．
	int i,				// in : [-]:		変換したい条番号．
	const Vector3d&x	// in : [m]:		変換したい慣性座標系．
) {
	Vector3d xyz = this->to_mycoord(x);
	Vector3d eta = this->SP[i].to_eta(xyz);
	return eta;
}

// 座標をナットの溝直交座標系から慣性座標系へ変換
Vector3d BS_Cylinder::to_inertialcoord
(						// out:	[m]:		慣性座標系の座標．
	int i,				// in : [-]:		変換したい条番号．
	const Vector3d&eta	// in :	[rad],[m]:	変換したい螺旋座標系．
) {
	Vector3d xyz = this->SP[i].to_xyz(eta);
	Vector3d x = this->to_inecoord(xyz);
	return x;
}

// ベクトルを慣性座標系からナットの溝直交座標系へ変換
Vector3d BS_Cylinder::to_etavelocity
(						// out:	[rad],[m]:	0成分：ナット位相角．1成分：eta座標．2成分：zeta座標．
	const Vector3d&v,	// in : [m/s]:		変換したい速度（螺旋座標系）．
	const Matrix3d&xyz2eta
) {
	Vector3d v_ = this->to_myvelocity(v);
	Vector3d v__ = xyz2eta * v_;
	return v__;
}

// ベクトルをナットの溝直交座標系から慣性座標系へ変換
Vector3d BS_Cylinder::to_inertialvelocity
(						// out:	[m]:		慣性座標系の座標．
	const Vector3d&v,	// in : [m/s]:		変換したい速度（螺旋座標系）．
	const Matrix3d&xyz2eta
) {
	Vector3d v_ = xyz2eta.transpose() * v;
	Vector3d v__ = this->to_inevelocity(v_);
	return v__;
}

// 螺旋位相角から変換行列をセットするメソッド．
Matrix3d BS_Cylinder::get_xyz2eta
(
	int i,				// in:	[-]:	螺旋番号
	double theta		// in:	[rad]:	現在の螺旋位相角
) {
	Matrix3d xyz2eta = this->SP[i].get_xyz2eta(theta);
	return xyz2eta;
}


// 慣性座標系のベクトルを螺旋座標系に変換するメソッド．
Vector3d BS_Cylinder::to_etavector(
	const Vector3d&a,	// in:	[any]	変換したいベクトル（螺旋座標系）
	const Matrix3d&xyz2eta
) {
	Vector3d b = this->to_myvector(a);
	Vector3d c = xyz2eta * b;

	return c;
}

// 螺旋座標系ベクトルを慣性座標系に変換するメソッド．
Vector3d BS_Cylinder::to_inertialvector(
	const Vector3d&a,	// in:	[any]	変換したいベクトル（螺旋座標系）
	const Matrix3d&xyz2eta
) {
	Matrix3d eta2xyz = xyz2eta.transpose();
	Vector3d b = eta2xyz * a;
	Vector3d c = this->to_inevector(b);
	return c;
}

// 接触楕円をトーラスに沿ってスライスするメソッド．
void BS_Cylinder::calc_slice
(
	int is,				// in:	[-]:	螺旋番号;
	int ig,				// in:	[-]:	溝番号．0 or 1;
	double th,			// in:	[rad]:	螺旋位相角．
	const Vector3d&p,	// in:	[m]:	接触点位置の中心．（慣性座標系）
	const Vector3d&xai,	// in:	[m]:	溝断面における螺旋進行方向ベクトル．（慣性座標系）
	double a,			// in:	[m]:	接触楕円長径．
	double b,			// in:	[m]:	接触楕円短径．
	const Vector2d*xy,	// in:	[-]:	スライス座標比 [-] xyは単位円の内側の値．0成分:a方向．1成分:b方向．
	int n,				// in:	[-]:	接触楕円のスライス数．
	const Vector2d&rho,	// in:	[1/m]:	各等価曲率．0成分:a方向．1成分:b方向．
	Vector3d*ps			// out:	[m]:	スライスされた各接触点位置．配列サイズ=n．（慣性座標系）
) {
	// 慣性座標系での，点pが乗るトーラスの垂直断面円の中心O．（トーラスをバナナに見立てると，バナナの芯のイメージ．）
	Vector3d O = this->to_inertialcoord(is,
		Vector3d(th, this->SP[is].GV[ig].eta[0], this->SP[is].GV[ig].eta[1]));

	// 円の計算はコストがかかるため，接触楕円が十分小さいことから2次近似を行う．
	Vector3d op = p - O;
	double op_norm = op.norm();

	Vector3d cn = op / op_norm; // 溝断面に垂直なベクトル
	Vector3d ct = cn.cross(xai);	// 接触点から見た溝Rの周方向．
	
	// 各スライスの中心点 pi の座標を算出
	for (int i = 0; i < n; i++) {
		double xi_ = a * xy[i].x();
		double yi_ = b * xy[i].y();
		double zi_ = - 0.5 * xi_ * xi_ * rho.x() - 0.5 * yi_ * yi_ * rho.y();
		ps[i] = p + cn * zi_ + ct * xi_ + xai * yi_;
	}
	return;
}

void BS_Cylinder::allocate(int ns) {

	this->nSP = ns;
	this->SP = new BS_Spiral[this->nSP];

	return;
};

// 入力ファイルをもとに総初期化を行う．
void BS_Cylinder::init(const BS_In::Cylinder & cylinder, const bool(&v_const)[3], const bool(&w_const)[3]) {

	this->set_const(v_const, w_const);

	this->rho = cylinder.density;
	this->nu = cylinder.poisson;
	this->E = cylinder.young;
	this->ri = cylinder.ri;
	this->ro = cylinder.ro;
	this->set_mI(cylinder.m, Vector3d(cylinder.Ix, cylinder.Iyz, cylinder.Iyz));

	for (int j = 0; j < this->nSP; j++)
		this->SP[j].init(cylinder.spiral[j]);

	return;
}

double BS_Cylinder::get_nd(int is) {
	return this->SP[is].get_nd();
}

double BS_Cylinder::get_r(int is) {
	return this->SP[is].get_r();
}

void BS_Cylinder::linspace
(
	int is,		// in:	[-]		条番号
	double th0,	// in:	[rad]	始点位相角		
	double th1,	// in:	[rad]	終点位相角
	int nb,		// in:	[-]		玉数
	Vector3d*xs	// out:	[m]		等配分された玉座標（慣性座標系）
) {
	VectorXd th_ = VectorXd::LinSpaced(nb, th0, th1);

	for (int i = 0; i < nb; i++)
		xs[i] = this->to_inertialcoord(is, Vector3d(th_[i], 0.0, 0.0));

	return;
}

void BS_Cylinder::save(BS_Out::Cylinder & OUT) {

	Vector3d ax = this->get_ax();
	for (int k = 0; k < 3; k++) {
		OUT.x[k] = this->x[k];
		OUT.v[k] = this->v[k];
		OUT.w[k] = this->w[k];
		OUT.ax[k] = ax[k];
		OUT.F[k] = this->F[k];
		OUT.T[k] = this->T[k];
		OUT.Fs[k] = this->Fs[k];
		OUT.Ts[k] = this->Ts[k];
	}
	Vector4d q = this->q.coeffs();
	for (int k = 0; k < 4; k++)
		OUT.q[k] = q[k];

	return;
}


BS_Cylinder::BS_Cylinder() {
	this->SP = NULL;
	return;
}

BS_Cylinder::~BS_Cylinder() {
	if (this->SP != NULL)
		delete[] this->SP;
	return;
}

