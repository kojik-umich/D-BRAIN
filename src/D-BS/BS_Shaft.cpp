/*******************************************************************************
!								"BS_Shaft.cpp"
!													2019/03/22	[Core-T]	楠崎
!	ねじ軸（以下シャフト）はモデリングの都合上，
!	1. 部材単体のシャフト（パーツ）
!	2. 単体シャフトを複数個組み合わせたもの（アセンブリ）
!	の二種類ある．（現状のアセンブリシャフトは1つの組み合わせしか考慮できない）
!	このクラスではアセンブリシャフトを定義しており，
!	Cylinderクラスでねじ溝の詳細な形状を定義している．
!
!*******************************************************************************/
#include "BS_Shaft.h"

void BS_Shaft::Nothing(void) {
}

void BS_Shaft::allocate(const std::vector<BS_In::Cylinder>&cylinders) {

	this->nCY = 1;
	this->CY = new BS_Cylinder[this->nCY];

	this->CY[0].allocate(cylinders[0].spiral.size());

	return;
}

void BS_Shaft::init(const std::vector<BS_In::Cylinder>&cylinders, const bool(&v_const)[3], const bool(&w_const)[3], const double(&ax0)[3], double v0, double w0) {

	this->CY[0].init(cylinders[0], v_const, w_const);

	this->rho = this->CY[0].rho = cylinders[0].density;
	this->nu = this->CY[0].nu = cylinders[0].poisson;
	this->E = this->CY[0].E = cylinders[0].young;

	this->w = this->get_ax() * w0;
	this->v = this->get_ax() * v0;

	Vector3d ax0_ = Vector3d(ax0) / ax0[0];		// x成分が1になるように一応規格化．
	this->tan_thy =-ax0_[2];
	this->tan_thz = ax0_[1];

	this->set_mI(cylinders[0].m, Vector3d(cylinders[0].Ix, cylinders[0].Iyz, cylinders[0].Iyz));

	this->set_const(v_const, w_const);

	//this->free_num = 6;
	//for (int i = 0; i < 3; i++) {
	//	this->free_num -= int(this->x_const[i]);
	//	this->free_num -= int(this->Rx_const[i]);
	//}
	return;
}

// 位置の初期化．とりあえず原点に置けばよい？
void BS_Shaft::init_pos(double v0, double w0) {

	this->x = Vector3d::Zero();

	Vector3d ax = Vector3d::UnitX();
	this->set_ax(ax);
	this->v = v0 * ax;
	this->w = w0 * ax;

	this->set_dx();

	return;
}

// step0（剛性計算）用のインタフェイス．
void BS_Shaft::get_y0(double*y0) {

	y0[0] = this->x.x() / Rigid::l;
	y0[1] = this->x.y() / Rigid::l;
	y0[2] = this->x.z() / Rigid::l;

	y0[3] = this->q.y();
	y0[4] = this->q.z();

	return;
}

// step0（剛性計算）用のインタフェイス．(拘束条件は要修正)
void BS_Shaft::set_y0(const double*y0, double v0, double w0) {

	Vector3d x_now = this->x;
	this->x = this->x_const.select(
		x_now,
		Vector3d(y0[0], y0[1], y0[2]) * Rigid::l
	);
	Vector3d q_ = this->Rx_const.select(
		Vector3d(this->q.x(), this->q.y(), this->q.z()),
		Vector3d(0.0, y0[3], y0[4])
	);
	this->q.y() = q_.y();
	this->q.z() = q_.z();
	this->q.normalize();

	Vector3d ax = this->get_ax();
	this->v = v0 * ax;
	this->w = w0 * ax;
	this->set_dx();

	return;
}

// シャフト（パーツ）の数値をシャフト（アセンブリ）の数値を基に更新．シャフト（アセンブリ）とシャフト（パーツ）は同一に動くものとする．
void BS_Shaft::set_dx(void) {
	this->CY[0].set_param(this->x, this->v, this->q, this->w);
	return;
}

// 値の設定メソッド．無次元量で渡された値を有次元にしてパラメタとして代入する．
void BS_Shaft::set_y_(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w) {

	this->set_y(x, v, q, w);
	this->set_dx();

	return;
}

// step0（剛性計算）用のインタフェイス．
void BS_Shaft::init_dyn0(void) {

	this->v = Vector3d::Zero();
	this->w = Vector3d::Zero();

	return;
}

void BS_Shaft::deinit_dyn0(double v0, double w0) {

	this->v = v0 * this->get_ax();
	this->w = w0 * this->get_ax();

	return;
}

// step0（剛性計算）用のインタフェイス．
void BS_Shaft::get_dyn_y0(double*y0) {

	for (int i = 0; i < 3; i++) {
		y0[i + 0] = this->x[i] / Rigid::l;
		y0[i + 3] = this->v[i] / Rigid::l * Rigid::t;
	}
	y0[6] = this->q.y() * 2;
	y0[7] = this->q.z() * 2;
	y0[8] = this->w.y() * Rigid::t;
	y0[9] = this->w.z() * Rigid::t;

	return;
}

// step0（剛性計算）用のインタフェイス．(拘束条件は要修正)
void BS_Shaft::set_dyn_y0(const double*y0) {

	Vector3d x_now = this->x;
	Vector3d x = this->x_const.select(x_now,
		Vector3d(y0[0], y0[1], y0[2]) * Rigid::l);
	Vector3d v = Vector3d(y0[3], y0[4], y0[5])
		* Rigid::l / Rigid::t;
	Vector3d w = Vector3d(0, y0[8], y0[9])
		/ Rigid::t;

	Vector3d t_now = 2 * Vector3d(0, this->q.y(), this->q.z());
	Vector3d t = this->Rx_const.select(t_now, Vector3d(0, y0[6], y0[7]));
	double   qy = 0.5 * t.y();
	double   qz = 0.5 * t.z();
	double   qw = sqrt(1.0 - qy * qy - qz * qz);

	this->x = x;
	this->v = v;
	this->q = Quaterniond(qw, 0, qy, qz);
	this->w = w;

	this->set_dx();

	return;
}

// F, T の入力により x, v, q, w 全ての時間微分値を返すメソッド．値を全て無次元量にして返す．
void BS_Shaft::get_dydt_(
	const Vector3d&F,	// in:	[N]:		外部荷重．
	const Vector3d&T,	// in:	[Nm]:		外部トルク．
	double dvdt0,		// in:	[m/s^2]:	速度の時間変化量
	double dwdt0,		// in:	[rad/s^2]	回転数の時間変化量
	double*dydt			// out:	[any]:		全微分値．
) {
	// 姿勢拘束（角速度とクォータニオンの両方を拘束，軸方向x回りの回転のみ物体座標系で拘束）
	Vector4d dqdt = this->get_dqdt_(this->Rx_const.y(), this->Rx_const.z(), this->tan_thy, this->tan_thz) * Rigid::t;
	Vector3d dwdt = this->get_dwdt_(T, this->Rx_const.y(), this->Rx_const.z(), dwdt0) * Rigid::t * Rigid::t;

	// 変位拘束(慣性座標系で拘束)
	Vector3d dxdt = this->get_dxdt() / Rigid::l * Rigid::t;
	Vector3d _dvdt = this->get_dvdt(F) / Rigid::l * Rigid::t * Rigid::t;
	Vector3d dvdt = (this->x_const == true).select(Vector3d::Zero(), _dvdt);
	if (this->x_const.x()) {
		dvdt += Vector3d(dvdt0, 0, 0);
	}

	dydt[0] = dxdt.x(); dydt[1] = dxdt.y(); dydt[2] = dxdt.z();
	dydt[3] = dvdt.x(); dydt[4] = dvdt.y(); dydt[5] = dvdt.z();
	dydt[6] = dqdt.w(); dydt[7] = dqdt.x(); dydt[8] = dqdt.y(); dydt[9] = dqdt.z();
	dydt[10] = dwdt.x(); dydt[11] = dwdt.y(); dydt[12] = dwdt.z();

	return;
}

// 角度拘束付きクォータニオンの時間微分式．y/z軸に回転拘束をかける場合，計算誤差で姿勢が徐々に変化するため，拘束をかけている．
Vector4d BS_Shaft::get_dqdt_(bool wy_const, bool wz_const, double tan_thy, double tan_thz) {
	// Quaterniond のコンストラクタは(w, x, y, z)の順番
	Quaterniond w_half(0.0, 0.5*this->w.x(), 0.5*this->w.y(), 0.5*this->w.z());
	Vector4d qw_half = (w_half * q).coeffs();	// クォータニオンは加算が定義されていないため，ベクトル配列で計算．
	Vector4d dqdt;
	double tau = 0.1;			// 緩和係数（本当は入力ファイルから入力できるようにするべき）
	double qx = q.x(), qy = q.y(), qz = q.z(), qw = q.w();

	// 完全拘束は角速度拘束だけで問題ないため，姿勢は拘束なし
	if (wy_const && wz_const) {
		dqdt = qw_half;
	}
	// wy方向が拘束されている場合，wy方向と軸方向に加速度成分を持たないようにする．
	// (= wy単位ベクトルと軸方向ベクトルの外積の向きのみに加速するようにする)
	else if (wy_const && !wz_const) {
		// Vector4d のコンストラクタは(x, y, z, w)の順番
		Vector4d _Cq(qz + tan_thy * qx, -qw - tan_thy * qy,
			qx - tan_thy * qz, -qy + tan_thy * qw);
		Vector4d Cq = _Cq * 2;
		double C = 2 * (qx * qz - qy * qw) +
			(qx * qx - qy * qy - qz * qz + qw * qw) * tan_thy;

		// 0割りの判定(多分いらない)
		double Cq_square = Cq.squaredNorm();
		if (Cq_square < 1e-20) {
			dqdt = qw_half;
		}
		else {
			dqdt = qw_half - Cq * (qw_half.dot(Cq) + 1 / tau * C) / Cq_square;
		}
	}
	// wz方向が拘束されている場合，wz方向と軸方向に加速度成分を持たないようにする．
	// (= wy単位ベクトルと軸方向ベクトルの外積の向きのみに加速するようにする)
	else if (!wy_const && wz_const) {
		// Vector4d のコンストラクタは(x, y, z, w)の順番
		Vector4d _Cq(qy - tan_thz * qx, qx + tan_thz * qy,
			qw + tan_thz * qz, qy - tan_thz * qw);
		Vector4d Cq = _Cq * 2;
		double C = 2 * (qx * qy + qz * qw) +
			(qx * qx - qy * qy - qz * qz + qw * qw) * tan_thz;

		// 0割りの判定（多分いらない）
		double Cq_sq = Cq.squaredNorm();
		if (Cq_sq < 1e-20) {
			dqdt = qw_half;
		}
		else {
			dqdt = qw_half - Cq * (qw_half.dot(Cq) + 1 / tau * C) / Cq_sq;
		}
	}
	// 軸方向の回転のみ拘束．
	else {
		dqdt = qw_half;
	}
	return dqdt;
}


// 角速度の時間微分式．軸自転方向と拘束されている成分は変化しないように設定．
Vector3d BS_Shaft::get_dwdt_(const Vector3d&T, bool wy_const, bool wz_const, double dwdt0) {
	Vector3d dwdt;

	// 完全拘束
	if (wy_const && wz_const) {
		Vector3d ax = this->get_ax();
		dwdt = ax * dwdt0;
	}
	// wy方向が拘束されている場合，wy方向と軸方向に加速度成分を持たないようにする．
	// (= wy単位ベクトルと軸方向ベクトルの外積の向きのみに加速するようにする)
	else if (wy_const && !wz_const) {
		Vector3d ax = this->get_ax();
		Vector3d ey = Vector3d(0, 1, 0);
		Vector3d ew = ax.cross(ey);			// y方向と軸方向に垂直なベクトル
		Vector3d _dwdt = this->get_dwdt(T);
		dwdt = _dwdt.dot(ew) * ew + ax * dwdt0;
	}
	// wz方向が拘束されている場合，wz方向と軸方向に加速度成分を持たないようにする．
	// (= wy単位ベクトルと軸方向ベクトルの外積の向きのみに加速するようにする)
	else if (!wy_const && wz_const) {
		Vector3d ax = this->get_ax();
		Vector3d ez = Vector3d(0, 0, 1);
		Vector3d ew = ax.cross(ez);			// z方向と軸方向に垂直なベクトル
		Vector3d _dwdt = this->get_dwdt(T);
		dwdt = _dwdt.dot(ew) * ew + ax * dwdt0;
	}
	// 軸方向の回転のみ拘束．
	else {
		Vector3d ax = this->get_ax();
		Vector3d dwdt_ = this->get_dwdt(T);
		dwdt = this->remove_ax(dwdt_) + ax * dwdt0;
	}
	return dwdt;
}


void BS_Shaft::set_Fall_Zero(void) {

	this->F = this->T = Vector3d::Zero();

	for (int i = 0; i < this->nCY; i++)
		this->CY[i].F = this->CY[i].T = this->CY[i].Fs = this->CY[i].Ts = Vector3d::Zero();

	return;
}

BS_Shaft::BS_Shaft() {
	this->CY = NULL;
	return;
}

BS_Shaft::~BS_Shaft() {
	if (this->CY != NULL)
		delete[] this->CY;
	return;
}





