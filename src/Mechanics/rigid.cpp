#include "Rigid.h"

double Rigid::l;
double Rigid::t;
Vector3d Rigid::g;

// コンストラクタ，軸方向ベクトル（定数値）を定義
Rigid::Rigid(void) {
	// 初期軸方向はx軸正を向くとする．
	this->ax0 << 1.0, 0.0, 0.0;
	this->ay0 << 0.0, 1.0, 0.0;
	this->az0 << 0.0, 0.0, 1.0;
	
	this->x_const = Eigen::Array<bool, 3, 1>(false, false, false);
	this->Rx_const = Eigen::Array<bool, 3, 1>(false, false, false);

	return;
}

// 値の設定メソッド．
void Rigid::set_param(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w) {
	this->x = x;
	this->v = v;
	this->q = q;
	this->w = w;
	this->q.normalize();
	return;
}

// 値の取得メソッド．
void Rigid::get_param(Vector3d&x, Vector3d&v, Quaterniond&q, Vector3d&w) {
	x = this->x;
	v = this->v;
	q = this->q;
	w = this->w;
	return;
}


// 値の設定メソッド．無次元量で渡された値を有次元にしてパラメタとして代入する．
void Rigid::set_y(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w) {
	this->x = x * Rigid::l;
	this->v = v * Rigid::l / Rigid::t;
	this->q = q;
	this->w = w / Rigid::t;
	this->q.normalize();
	return;
}


// 値の取得メソッド．有次元量で保持している値を無次元量にして返すメソッド．
void Rigid::get_y(Vector3d&x, Vector3d&v, Quaterniond&q, Vector3d&w) {
	x = this->x / Rigid::l;
	v = this->v / Rigid::l * Rigid::t;
	q = this->q;
	w = this->w * Rigid::t;
	return;
}

// そのまんまだけど現在の速度を出力する．
Vector3d Rigid::get_dxdt() {
	return this->v;
}

// 質量から微小速度変化量を出力する．
Vector3d Rigid::get_dvdt(const Vector3d&F) {
	Vector3d dvdt = this->m_inv * F;
	return dvdt;
}

// 現在のクォータニオンから角速度の入力により微小クォータニオン変化量を出力する．
// クォータニオン同士の外積として計算
Quaterniond Rigid::get_dqdt(void) {
	Quaterniond w_half(0.0, 0.5*this->w.x(), 0.5*this->w.y(), 0.5*this->w.z());
	Quaterniond dqdt = w_half * this->q;
	return dqdt;
}

// 現在の角速度からトルクの入力により微小角速度変化量を出力する．
Vector3d Rigid::get_dwdt(const Vector3d&T) {
	Vector3d Iw = this->I.cwiseProduct(this->w);
	Vector3d wIw = this->w.cross(Iw);
	Vector3d dwdt = this->I_inv.cwiseProduct(T - wIw);

	return dwdt;
}

// F, T の入力により x, v, q, w 全ての時間微分値を返すメソッド．拘束条件付き．値を全て無次元量にして返す．
void Rigid::get_dydt(
	const Vector3d&F,	// in:	[N]:	外部荷重．
	const Vector3d&T,	// in:	[Nm]:	外部トルク．
	double*dydt			// out:	[-]:	全微分値（無次元量）．
) {
	Vector3d	dxdt = this->get_dxdt() / Rigid::l * Rigid::t;
	Quaterniond	dqdt = this->get_dqdt();
	Vector3d dvdt = this->get_dvdt(F) / Rigid::l * Rigid::t * Rigid::t; 
	Vector3d dwdt = this->get_dwdt(T) * Rigid::t * Rigid::t;

	dvdt = (this->x_const == true).select(Vector3d::Zero(), dvdt); // ※クォータニオンに拘束条件をかけていないため，正しく拘束されない
	dwdt = (this->Rx_const == true).select(Vector3d::Zero(), dwdt);


	dydt[0]  = dxdt.x(); dydt[1]  = dxdt.y(); dydt[2]  = dxdt.z();
	dydt[3]  = dvdt.x(); dydt[4]  = dvdt.y(); dydt[5]  = dvdt.z();
	dydt[6]  = dqdt.w() * Rigid::t; dydt[7]  = dqdt.x() * Rigid::t; dydt[8]  = dqdt.y() * Rigid::t; dydt[9]  = dqdt.z() * Rigid::t;
	dydt[10] = dwdt.x(); dydt[11] = dwdt.y(); dydt[12] = dwdt.z();
	return;
}

// 拘束条件定義
void Rigid::set_const(const bool (&v_const)[3], const bool (&w_const)[3]) {
	for (int i = 0; i < 3; i++) {
		this->x_const[i] = v_const[i];
		this->Rx_const[i] = w_const[i];
	}
	return;
}

// 入力された慣性系座標を，自分の系の座標に変換するメソッド．
Vector3d Rigid::to_mycoord(const Vector3d&x) {
	Quaterniond qinv = this->q.inverse();
	Vector3d dx = x - this->x;
	Vector3d y  = qinv * dx;
	return y;
}

// 入力された自分の系の座標を，慣性系座標に変換するメソッド．
Vector3d Rigid::to_inecoord(const Vector3d&y) {
	Vector3d dx = this->q * y;
	Vector3d x  = dx + this->x;
	return x;
}

// 入力された慣性系速度を，自分の系の速度に変換するメソッド．
Vector3d Rigid::to_myvelocity(const Vector3d&v) {
	Quaterniond qinv = this->q.inverse();
	Vector3d dv = v - this->v;
	Vector3d V  = qinv * dv;
	return V;
}

// 入力された自分の系の座標を，慣性系座標に変換するメソッド．
Vector3d Rigid::to_inevelocity(const Vector3d&V) {
	Vector3d dv = this->q * V;
	Vector3d v  = dv + this->v;
	return v;
}

// 入力された慣性系ベクトルを自分の系のベクトルに変換するメソッド．
Vector3d Rigid::to_myvector(const Vector3d&x) {
	Quaterniond qinv = this->q.inverse();
	Vector3d y  = qinv * x;
	return y;
}

// 入力された自分の系のベクトルを慣性系ベクトルに変換するメソッド．
Vector3d Rigid::to_inevector(const Vector3d&y) {
	Vector3d x  = this->q * y;
	return x;
}

// 入力された慣性座標系上の点について，その速度を求める（自分の体内の点かは判定しない）
Vector3d Rigid::surface_velocity(const Vector3d&x) {
	Vector3d r = x - this->x;
	Vector3d V = w.cross(r);
	V = V + this->v;
	return V;
}


// 入力された慣性座標系上の点にかかる力でトルクを求める．（自分の体内の点かは判定しない）
Vector3d Rigid::calc_Torque(const Vector3d&x, const Vector3d&F) {
	Vector3d r = x - this->x;
	Vector3d T = r.cross(F);
	return T;
}


// 入力された位置・方向に力がかかった際のトルクの向きを求める．（自分の体内の点かは判定しない）
Vector3d Rigid::calc_TorqueDirection(const Vector3d&x, const Vector3d&u) {
	Vector3d r = x - this->x;
	Vector3d T = r.cross(u);
	return T.normalized();
}


// 現在の軸方向を算出するメソッド．（慣性座標系，ノルム=1）
Vector3d Rigid::get_ax(void) {
	Vector3d ax = this->to_inevector(this->ax0);
	return ax;
}


// 現在の方位角180°の径方向を算出するメソッド．造語．（慣性座標系，ノルム=1）
Vector3d Rigid::get_ay(void) {
	Vector3d ay = this->to_inevector(this->ay0);
	return ay;
}


// 現在の方位角270°の径方向を算出するメソッド．造語．（慣性座標系，ノルム=1）
Vector3d Rigid::get_az(void) {
	Vector3d az = this->to_inevector(this->az0);
	return az;
}


// 入力された重力加速度ベクトルから自身の質量の積をとり，荷重ベクトルを算出する．
Vector3d Rigid::get_mg(void) {
	return this->m * Rigid::g;
}

// 入力されたベクトルの，軸方向成分を除去し，軸ベクトルと直角なベクトルに成形する．
Vector3d Rigid::remove_ax(const Vector3d&x) {
	Vector3d ax = this->get_ax();
	return x - ax.dot(x) * ax;
}

//　入力した軸方向を向くよう，Quaternionを計算し設定する．q.x()は0で固定．
void Rigid::set_ax(const Vector3d&ax) {
	Vector3d ax_ = ax.normalized();
	double y2z2 = 0.5 * (1.0 - ax_.x());
	double sq_1_y2z2 = sqrt(1 - y2z2);
	double y = -0.5 * ax_.z() / sq_1_y2z2;
	double z =  0.5 * ax_.y() / sq_1_y2z2;

	Quaterniond q_= Quaterniond(sq_1_y2z2, 0.0, y, z);
	this->q = q_;

}

// トルクと荷重をメンバ変数に代入
void Rigid::set_FT(const Vector3d & F, const Vector3d & T) {
	this->F = F;
	this->T = T;
	return;
}

void Rigid::set_mI(double m, const Vector3d & I) {
	this->m = m;
	this->m_inv = 1.0 / m;
	this->I = I;
	this->I_inv = I.cwiseInverse();
	return;
}

// x,v,q,w の4パラメータを外に書き出すメソッド．（13変数）
void Rigid::save(double*x, double*v, double*q, double*w, double*ax, double*F, double*T) {

	Vector3d ax_ = this->get_ax();
	for (int i = 0; i < 3; i++) {
		x[i] = this->x[i];
		v[i] = this->v[i];
		w[i] = this->w[i];
		ax[i] = ax_[i];
		F[i] = this->F[i];
		T[i] = this->T[i];
	}
	q[0] = this->q.x();
	q[1] = this->q.y();
	q[2] = this->q.z();
	q[3] = this->q.w();
	return;
}

Rigid::~Rigid(void) {
}
