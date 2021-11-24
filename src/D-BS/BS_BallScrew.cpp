/*******************************************************************************
!								"BS_BallScrew.cpp"
!													2020/04/15	[Core-T]	楠崎
!
! ボールねじオブジェクト．他で言うところの軸受オブジェクトに相当．
!
!*******************************************************************************/

#include "BS_BallScrew.h"

// ボールねじ構成要素を動的確保する．一番最初に起動してほしいメソッド．
void BS_BallScrew::allocate(const BS_In&IN) {

	if (IN.nut.size() == 1)
		this->NT = new BS_SingleNut();
	else
		this->NT = new BS_DoubleNut();

	this->nCC = IN.circuit.size();
	this->CC = new BS_Circuit[this->nCC];
	for (size_t i = 0; i < IN.circuit.size(); i++)
		this->CC[i].allocate(IN.circuit[i].ball.size());

	this->nP = 0;
	for (size_t i = 0; i < IN.circuit.size(); i++)
		this->nP += IN.circuit[i].ball.size();

	this->BNP = new BS_BallNutPair[this->nP];
	this->BSP = new BS_BallShaftPair[this->nP];

	this->NT->allocate(IN.nut);
	this->ST.allocate(IN.shaft);

	return;
}

// ボールと対応する部材の番号を紐付ける．
void BS_BallScrew::link(const BS_In&IN) {

	int iter = 0;

	for (int i = 0; i < this->nCC; i++) {

		this->CC[i].link(
			&this->NT->CY[IN.circuit[i].inut],
			IN.circuit[i].is);

		for (int j = 0; j < this->CC[i].nBL; j++) {

			this->BNP[iter].link(
				&this->CC[i].BL[j],
				&this->NT->CY[IN.circuit[i].inut],
				IN.circuit[i].is);


			this->BSP[iter].link(
				&this->CC[i].BL[j],
				&this->ST.CY[0],
				IN.circuit[i].is);
			iter++;
		}
	}
	return;
}

// 回路にある各玉の座標を初期化するメソッド．とりあえず0始まりで試してみる．
void BS_BallScrew::init_position(
	double wn,	// in:	[rad/s]		理論公転数．
	double ws	// in:	[rad/s]		理論公転数．
) {
	for (int i = 0; i < this->nCC; i++) {

		Vector3d dx = this->CC[i].CY->x;

		Vector3d*x0s = new Vector3d[this->CC[i].nBL];
		this->CC[i].get_x0s(x0s);

		for (int j = 0; j < this->CC[i].nBL; j++)
			this->CC[i].BL[j].x = x0s[j] + 0.5 * dx;

		double nd = this->CC[i].get_nd();
		double pcd = this->CC[i].get_r() * 2;

		for (int j = 0; j < this->CC[i].nBL; j++) {
			Vector3d t = this->CC->get_t(j);
			double D = this->CC[i].BL[j].r * sqrt(2.0);	// 接触角45度を仮定．
			double w = 0.5 * (1.0 - D / pcd) * ws + 0.5 * (1.0 + D / pcd) * wn;	// 公転角速度
			Vector3d v = w * nd * t;
			Vector3d w_ = (pcd / D - D / pcd) * 0.5 * (wn - ws) * Vector3d(1.0, 0.0, 0.0);	// 自転角速度

			this->CC[i].BL[j].v = v;
			this->CC[i].BL[j].w = w_;
			this->CC[i].BL[j].q.setIdentity();
		}
	}
	return;
}
// 荷重の初期化．また，荷重からモーメントを計算し，入力モーメントに加算．
void BS_BallScrew::init_Load(const BS_In&IN) {

	Ball Dummy;
	Dummy.x = Vector3d::Zero();
	this->LD.F = Vector3d::Zero();
	this->LD.T = Vector3d::Zero();

	for (size_t i = 0; i < IN.load.size(); i++) {

		this->LD.F += Vector3d(IN.load[i].F);
		this->LD.T += Dummy.calc_Torque(
			Vector3d(IN.load[i].x),
			Vector3d(IN.load[i].F))
			+ Vector3d(IN.load[i].T);
	}
	return;
}

// 入力からボールねじの初期化を行う．
void BS_BallScrew::init(const BS_In&IN, double v0, double w0, double wn) {

	this->link(IN);

	for (int i = 0; i < this->nCC; i++)
		this->CC[i].init(IN.circuit[i]);

	this->NT->init(IN.nut, wn);
	this->ST.init(IN.shaft, IN.bound.v_const, IN.bound.w_const, IN.bound.tan_thy, IN.bound.tan_thz, v0, w0);

	for (int i = 0; i < this->nP; i++) {
		this->BNP[i].init(IN.BallNutPair[i], IN.tribology, IN.oil);
		this->BSP[i].init(IN.BallShaftPair[i], IN.tribology, IN.oil);
	}
	this->init_position(wn, w0);
	this->ST.init_pos(v0, w0);

	for (int i = 0; i < this->nP; i++) {
		Vector3d etan = this->BNP[i].get_eta();
		this->BNP[i].th0 = etan[0];
		Vector3d etas = this->BSP[i].get_eta();
		this->BSP[i].thp = etas[0];
	}
	this->init_Load(IN);

	for (int j = 0; j < this->nCC; j++) {
		this->CC[j].iSP = IN.circuit[j].is;
		this->CC[j].nBL = IN.circuit[j].ball.size();

		this->CC[j].th0 = IN.circuit[j].th0;
		this->CC[j].th1 = IN.circuit[j].th1;
	}
	return;
}

// 【初期値計算】シャフト，ボールの初期位置をざっくり決めるメソッド．
void BS_BallScrew::preset_y0(double dx0, double dth0, double dx1) {

	this->preset_y0_F(dx0);
	this->preset_y0_T(dth0);
	this->preset_y0_x(dx1);

	return;
}

// 入力にある初期位置にシャフトを設定するメソッド．
void BS_BallScrew::lock_y0(const double*x0, const double*ax0, double v0, double w0) {

	Vector3d x = Vector3d(x0);
	Vector3d ax = Vector3d(ax0);

	this->ST.x = x;
	this->ST.set_ax(ax);
	this->ST.v = v0 * ax;
	this->ST.w = w0 * ax;

	this->ST.set_dx();

	for (int i = 0; i < this->nCC; i++)
		for (int j = 0; j < this->CC[i].nBL; j++)
			this->CC[i].BL[j].x += Vector3d(x0) / 2;

	return;
}

// 【初期値計算】ボール全体とシャフトを荷重方向に動かして，ざっくり初期位置を決めるメソッド．
void BS_BallScrew::preset_y0_F(
	double dx0		// in:	[m]	初期微小変化量．通常1nmなど系に対して十分小さい値を用いてください．
) {
	// 荷重が十分小さい場合，このメソッドは使われない．
	double F_norm = this->LD.F.norm();
	if (F_norm < 1e-20)
		return;

	Vector3d e = this->LD.F.normalized();
	double dx, F0, F1;
	F0 = 0.0;
	dx = dx0;

	// まずボール全体を動かす．（移動量は指数関数増加）
	for (int i = 0; i < 1e3; i++) {
		// ボール全体を荷重方向に動かす．
		for (int i = 0; i < this->nCC; i++)
			for (int j = 0; j < this->CC[i].nBL; j++)
				this->CC[i].BL[j].x += e * dx;

		// ナットにかかる荷重を求める．
		Vector3d Fn_sum = Vector3d::Zero();
		for (int j = 0; j < this->nP; j++) {
			Vector2d Fbn;
			Vector3d Fnb, Tnb;
			this->BNP[j].get_F0(Fbn, Fnb, Tnb);
			Fn_sum += Fnb;
		}
		// 荷重を更新．
		F1 = Fn_sum.norm();

		// 荷重が外部荷重を超過している場合，ループを抜ける．
		if (F1 > F_norm) {
			// 線形補間でおおよその位置に戻す．

			for (int i = 0; i < this->nCC; i++)
				for (int j = 0; j < this->CC[i].nBL; j++)
					this->CC[i].BL[j].x -= (F1 - F_norm) / (F1 - F0) * e * dx;
			break;
		}
		// 補正値を更新．
		dx *= 1.1;
		F0 = F1;
	}
	// 次にシャフトを動かす．（移動量は指数関数増加）
	F0 = 0.0;
	dx *= 0.1;
	for (int i = 0; i < 1e3; i++) {
		// シャフトを荷重方向に動かす．
		this->ST.x += e * dx;
		this->ST.set_dx();

		// シャフトにかかる荷重を求める．
		Vector3d Fs_sum = Vector3d::Zero();
		for (int j = 0; j < this->nP; j++) {
			Vector2d Fbs;
			Vector3d Fsb, Tsb;
			this->BSP[j].get_F0(Fbs, Fsb, Tsb);
			Fs_sum += Fsb;
		}
		// 荷重を更新．
		F1 = Fs_sum.norm();

		// 荷重が外部荷重を超過している場合，ループを抜ける．
		if (F1 > F_norm) {
			// 線形補間でおおよその位置に戻す．
			this->ST.x -= (F1 - F_norm) / (F1 - F0) * e * dx;
			this->ST.set_dx();
			break;
		}
		// 補正値を更新．
		dx *= 1.1;
		F0 = F1;
	}
	return;
}

// シャフトをモーメント方向に回転して，玉と接触するように初期位置を決定
void BS_BallScrew::preset_y0_T(double dth0) {

	// 荷重が十分小さい場合，このメソッドは使われない．
	Vector2d T2d(this->LD.T.y(), this->LD.T.z());
	double T_norm = 2 * T2d.norm();
	if (T_norm < 1e-20)
		return;

	Vector2d e = T2d.normalized();
	double dth, T0, T1;
	T0 = 0.0;
	dth = dth0;

	// シャフトを移動（移動量は指数関数的に増加）
	for (int i = 0; i < 1e3; i++) {

		Vector3d ax = this->ST.get_ax();
		this->ST.set_ax(Vector3d(1.0, ax.y() + dth * e[1], ax.z() - dth * e[0]));
		this->ST.set_dx();

		Vector2d Ts_sum = Vector2d::Zero();
		for (int j = 0; j < this->nP; j++) {
			Vector2d Fbs;
			Vector3d Fsb, Tsb;
			this->BSP[j].get_F0(Fbs, Fsb, Tsb);
			Ts_sum += Vector2d(Tsb.y(), Tsb.z());
		}
		// 荷重を更新．
		T1 = Ts_sum.norm();

		// 荷重が外部荷重を超過している場合，ループを抜ける．
		if (T1 > T_norm) {
			// 線形補間でおおよその位置に戻す．
			double dth_ = -(T1 - T_norm) / (T1 - T0) * dth;
			Vector3d ax = this->ST.get_ax();
			this->ST.set_ax(Vector3d(1.0, ax.y() + dth_ * e[1], ax.z() - dth_ * e[0]));
			this->ST.set_dx();
			break;
		}
		// 補正値を更新．
		dth *= 1.1;
		T0 = T1;
	}
	return;
}

// ナット側にボールが触れるまでくっつけるメソッド．
void BS_BallScrew::preset_y0_x(double dx) {

	for (int i = 0; i < this->nP; i++) {
		// まず螺旋のもとの位置に戻してやる．
		Vector3d eta_zeta = this->BNP[i].get_eta();
		double eta = eta_zeta[1];
		double zeta = eta_zeta[2];
		this->BNP[i].set_eta0(Vector2d(eta, zeta));

		for (int j = 0; j < 1e3; j++) {
			Vector2d Fbn;
			Vector3d Fnb, Tnb;
			this->BNP[i].get_F0(Fbn, Fnb, Tnb);
			double F = abs(Fbn.x());

			// 接していたら抜ける．
			if (F > 1e-20)
				break;

			eta -= dx;
			dx *= 1.1;
			this->BNP[i].set_eta0(Vector2d(eta, zeta));
		}
	}
	return;
}

// step0（剛性計算）用のインタフェイス．
void BS_BallScrew::get_y0(double*y0) {

	this->ST.get_y0(y0);

	for (int i = 0; i < this->nP; i++) {
		int i2 = i * 2;
		Vector3d eta = this->BNP[i].get_eta();
		y0[i2 + 5] = eta[1] / Rigid::l;
		y0[i2 + 6] = eta[2] / Rigid::l;
	}
	return;
}

// step0（剛性計算）：変数y0を入力
void BS_BallScrew::set_y0(
	const double*y0,	// in :[-]		: 変数y0
	double v0,			// in :[m/s]	: シャフト進行速度
	double w0			// in :[rad/s]	: シャフト回転速度
) {

	this->ST.set_y0(y0, v0, w0);

	for (int i = 0; i < this->nP; i++) {
		int i2 = i * 2;
		Vector2d eta = Vector2d(
			y0[i2 + 5] * Rigid::l,
			y0[i2 + 6] * Rigid::l
		);
		this->BNP[i].set_eta0(eta);
	}
	return;
}

// 荷重を配列にして返すメソッド．
void BS_BallScrew::get_F0(double*F0) {

	Vector3d Fs = this->LD.F;
	Vector3d Ts = this->LD.T;

	for (int ib = 0; ib < this->nP; ib++) {

		Vector2d Fbn, Fbs;
		Vector3d Fnb, Tnb, Fsb, Tsb;

		this->BNP[ib].get_F0(Fbn, Fnb, Tnb);
		this->BSP[ib].get_F0(Fbs, Fsb, Tsb);

		Vector2d Fb = Fbn + Fbs;

		int i2 = ib * 2;

		// 弱いバネ（斥力）．x^-1 の形にすることによって，外側に収束しやすくする．
		double x = Vector2d(this->BNP[ib].BL->x[1], this->BNP[ib].BL->x[2]).norm();
		double k = 1e-20;

		F0[i2 + 5] = Fb[0] - k / x;
		F0[i2 + 6] = Fb[1];

		Fs += Fsb;
		Ts += Tsb;
	}
	Vector3d Ts_ = Ts / Rigid::l;

	return;
}

// 荷重を配列にして返すメソッド．
void BS_BallScrew::get_F1(double*F1) {

	Vector3d Fs = this->LD.F;
	Vector3d Ts = this->LD.T;

	bool is_RightHand = (this->NT->w.x() - this->ST.w.x()) > 0;	// シャフトから見て玉が反時計回りに公転しているときtrue

	for (int ib = 0; ib < this->nP; ib++) {

		Vector2d Fbn, Fbs;
		Vector3d Fnb, Fsb, Tnb, Tsb;
		this->BNP[ib].get_F1(!is_RightHand, Fbn, Fnb, Tnb);
		this->BSP[ib].get_F1(is_RightHand, Fbs, Fsb, Tsb);
		Vector2d Fb = Fbn + Fbs;
		//double Fbn_ = Fbn.norm();
		//double Fbs_ = Fbs.norm();

		//double sin_th = (Fbn[1] * Fbs[0] - Fbn[0] * Fbs[1]) / (Fbn_ * Fbs_);

		int i2 = ib * 2;
		F1[i2 + 5] = Fb[0]; // Fbn_ - Fbs_;
		F1[i2 + 6] = Fb[1]; // sin_th * sqrt(Fbn_ * Fbs_);

		Fs += Fsb;
		Ts += Tsb;
	}
	Vector3d Ts_ = Ts / Rigid::l;

	F1[0] = Fs[0];
	F1[1] = Fs[1];
	F1[2] = Fs[2];
	F1[3] = Ts_[1];
	F1[4] = Ts_[2];

	return;
}

// 各転動体の純転がり速度を求めるメソッド．
void BS_BallScrew::pure_Rolling(void) {

	// 全ての転動体を純転がり状態に設定する．
	for (int i = 0; i < this->nP; i++) {

		// 各溝での接触点数を求める．
		int nn = this->BNP[i].num_Contact();
		int ns = this->BSP[i].num_Contact();
		std::cout << i << std::endl;
		std::cout << nn << std::endl;
		std::cout << ns << std::endl << std::endl;
	}
}

void BS_BallScrew::get_y2(double * y2) {

	for (int ib = 0; ib < this->nP; ib++) {

		Vector3d v = this->BNP[ib].BL->v / Rigid::t * Rigid::l;
		Vector3d w = this->BNP[ib].BL->w / Rigid::t;

		for (int j = 0; j < 3; j++) {
			y2[j + 0] = v[j];
			y2[j + 3] = w[j];
		}
	}
	return;
}

void BS_BallScrew::set_y2(const double * y2) {

	for (int ib = 0; ib < this->nP; ib++) {

		Vector3d v(y2[0], y2[1], y2[2]);
		Vector3d w(y2[3], y2[4], y2[5]);
		this->BNP[ib].BL->v = v * Rigid::t / Rigid::l;
		this->BNP[ib].BL->w = w * Rigid::t;
	}
	return;
}

void BS_BallScrew::get_F2(double * f2) {

	for (int ib = 0; ib < this->nP; ib++) {

		Vector3d vFn, vTn;
		this->BNP[ib].get_F2(vFn, vTn);

		Vector3d vFs, vTs;
		this->BSP[ib].get_F2(vFs, vTs);

		Vector3d vF = vFn + vFs;
		Vector3d vT = vTn + vTs;

		for (int j = 0; j < 3; j++) {
			f2[j + 0] = vF[j];
			f2[j + 3] = vT[j];
		}
	}
	return;
}


void BS_BallScrew::init_dyn0(void) {

	this->ST.init_dyn0();

	for (int ib = 0; ib < this->nP; ib++) {
		this->BNP[ib].mem_BLv = this->BNP[ib].BL->v;
		this->BNP[ib].BL->v.setZero();
	}
	return;
}

void BS_BallScrew::deinit_dyn0(void) {

	for (int ib = 0; ib < this->nP; ib++)
		this->BNP[ib].BL->v = this->BNP[ib].mem_BLv;

	this->ST.deinit_dyn0();

	return;
}

void BS_BallScrew::get_dyn_y0(double * y0) {

	this->ST.get_dyn_y0(y0);

	for (int ib = 0; ib < this->nP; ib++) {
		int i5 = ib * 5;
		Vector3d eta = this->BNP[ib].get_eta0();
		y0[i5 + 11] = eta[1] / Rigid::l;
		y0[i5 + 12] = eta[2] / Rigid::l;
		for (int j = 0; j < 3; j++)
			y0[i5 + 13 + j] = this->BNP[ib].BL->v[j]
			/ Rigid::l * Rigid::t;
	}
	return;
}

// step0（剛性計算）：変数y0を入力
void BS_BallScrew::set_dyn_y0(const double*y0) {

	this->ST.set_dyn_y0(y0);

	for (int ib = 0; ib < this->nP; ib++) {
		int i5 = ib * 5;
		Vector2d eta = Vector2d(
			y0[i5 + 11] * Rigid::l,
			y0[i5 + 12] * Rigid::l
		);
		this->BNP[ib].set_eta0(eta);
		this->BNP[ib].BL->v =
			Vector3d(y0[i5 + 13], y0[i5 + 14], y0[i5 + 15])
			* Rigid::l / Rigid::t;
	}
	return;
}

// 荷重を配列にして返すメソッド．
void BS_BallScrew::get_dyn_dydt0(double*dydt) {

	Vector3d Fst = this->LD.F;
	Vector3d Tst = this->LD.T;

	bool is_RightHand = (this->NT->w.x() - this->ST.w.x()) > 0;	// シャフトから見て玉が反時計回りに公転しているときtrue

	for (int ib = 0; ib < this->nP; ib++) {

		int i5 = ib * 5;

		Vector3d etavn, etavs, Fbn, Fbs, Fnb, Fsb, Tnb, Tsb;
		this->BNP[ib].get_dyn_F0(!is_RightHand, etavn, Fbn, Fnb, Tnb);
		this->BSP[ib].get_dyn_F0(is_RightHand, etavs, Fbs, Fsb, Tsb);
		Vector3d Fb = Fbn + Fbs;

		for (int j = 0; j < 2; j++)
			dydt[i5 + 11 + j] = etavn[j + 1];

		for (int j = 0; j < 3; j++)
			dydt[i5 + 13 + j] = Fb[j] * this->BNP[ib].BL->m_inv;

		Fst += Fsb;
		Tst += Tsb;
	}

	for (int i = 0; i < 3; i++) {
		dydt[i + 0] = this->ST.v[i];
		dydt[i + 3] = Fst[i] * this->ST.m_inv;
		dydt[i + 8] = Tst[i] * this->ST.I_inv[i];
	}
	for (int i = 0; i < 2; i++)
		dydt[i + 6] = this->ST.w[i + 1];

	return;
}

void BS_BallScrew::set_dyn_y1(const double*y) {

	Vector3d x, v, w; Quaterniond q;
	x.x() = y[0];
	x.y() = y[1];
	x.z() = y[2];
	v.x() = y[3];
	v.y() = y[4];
	v.z() = y[5];
	q.w() = y[6];
	q.x() = y[7];
	q.y() = y[8];
	q.z() = y[9];
	w.x() = y[10];
	w.y() = y[11];
	w.z() = y[12];
	this->NT->set_y_(x, v, q, w);

	x.x() = y[13];
	x.y() = y[14];
	x.z() = y[15];
	v.x() = y[16];
	v.y() = y[17];
	v.z() = y[18];
	q.w() = y[19];
	q.x() = y[20];
	q.y() = y[21];
	q.z() = y[22];
	w.x() = y[23];
	w.y() = y[24];
	w.z() = y[25];
	this->ST.set_y_(x, v, q, w);

	for (int i = 0; i < this->nP; i++) {
		int i13 = i * 13 + 26;
		x.x() = y[i13 + 0];
		x.y() = y[i13 + 1];
		x.z() = y[i13 + 2];
		v.x() = y[i13 + 3];
		v.y() = y[i13 + 4];
		v.z() = y[i13 + 5];
		q.w() = y[i13 + 6];
		q.x() = y[i13 + 7];
		q.y() = y[i13 + 8];
		q.z() = y[i13 + 9];
		w.x() = y[i13 + 10];
		w.y() = y[i13 + 11];
		w.z() = y[i13 + 12];
		this->BNP[i].BL->set_y(x, v, q, w);
	}
	return;
}

void BS_BallScrew::get_dyn_y1(double*y) {

	Vector3d x, v, w; Quaterniond q;
	this->NT->get_y(x, v, q, w);
	y[0] = x.x();
	y[1] = x.y();
	y[2] = x.z();
	y[3] = v.x();
	y[4] = v.y();
	y[5] = v.z();
	y[6] = q.w();
	y[7] = q.x();
	y[8] = q.y();
	y[9] = q.z();
	y[10] = w.x();
	y[11] = w.y();
	y[12] = w.z();

	this->ST.get_y(x, v, q, w);
	y[13] = x.x();
	y[14] = x.y();
	y[15] = x.z();
	y[16] = v.x();
	y[17] = v.y();
	y[18] = v.z();
	y[19] = q.w();
	y[20] = q.x();
	y[21] = q.y();
	y[22] = q.z();
	y[23] = w.x();
	y[24] = w.y();
	y[25] = w.z();

	for (int i = 0; i < this->nP; i++) {
		int i13 = i * 13 + 26;
		this->BNP[i].BL->get_y(x, v, q, w);
		y[i13 + 0] = x.x();
		y[i13 + 1] = x.y();
		y[i13 + 2] = x.z();
		y[i13 + 3] = v.x();
		y[i13 + 4] = v.y();
		y[i13 + 5] = v.z();
		y[i13 + 6] = q.w();
		y[i13 + 7] = q.x();
		y[i13 + 8] = q.y();
		y[i13 + 9] = q.z();
		y[i13 + 10] = w.x();
		y[i13 + 11] = w.y();
		y[i13 + 12] = w.z();
	}
	return;
}

void BS_BallScrew::get_dyn_dydt1(
	double*dydt,	// out:	微分配列
	double dvdt,	// in:	
	double dwdt 	// in:	
) {
	// 各玉の接触力を求める．
	this->NT->F = this->ST.F = this->NT->T = this->ST.T = Vector3d::Zero();
	for (int i = 0; i < this->NT->nCY; i++)
		this->NT->CY[i].F = this->NT->CY[i].T = this->NT->CY[i].Fs = this->NT->CY[i].Ts = Vector3d::Zero();
	for (int i = 0; i < this->ST.nCY; i++)
		this->ST.CY[i].F = this->ST.CY[i].T = this->ST.CY[i].Fs = this->ST.CY[i].Ts = Vector3d::Zero();

	for (int i = 0; i < this->nP; i++) {

		Vector3d Fbn, Tbn, Tnb, Fnbs, Tnbs;
		this->BNP[i].get_FT(Fbn, Tbn, Tnb, Fnbs, Tnbs);
		this->BNP[i].CY->F += -Fbn;
		this->BNP[i].CY->T += Tnb;
		this->BNP[i].CY->Fs += Fnbs;
		this->BNP[i].CY->Ts += Tnbs;

		Vector3d Fbs, Tbs, Tsb, Fsbs, Tsbs;
		this->BSP[i].get_FT(Fbs, Tbs, Tsb, Fsbs, Tsbs);
		this->BSP[i].CY->F += -Fbs;
		this->BSP[i].CY->T += Tsb;
		this->BSP[i].CY->Fs += Fsbs;
		this->BSP[i].CY->Ts += Tsbs;

		// 重力加速度を加算．
		this->BNP[i].BL->F = Fbn + Fbs + this->BNP[i].BL->get_mg();
		this->BNP[i].BL->T = Tbn + Tbs;
		double dydt_[13];
		this->BNP[i].BL->get_dydt(this->BNP[i].BL->F, this->BNP[i].BL->T, dydt_);

		// ループごとに戻り値に代入．
		int i13 = i * 13 + 26;
		for (int j = 0; j < 13; j++)
			dydt[i13 + j] = dydt_[j];

		this->NT->F += -Fbn;
		this->NT->T += Tnb;
		this->ST.F += -Fbs;
		this->ST.T += Tsb;
	}
	// 一応微分をとるが，NTはv_const=trueなので0が戻り値となる．
	double dydtn[13];
	this->NT->get_dydt(this->NT->F, this->NT->T, dydtn);
	for (int j = 0; j < 13; j++)
		dydt[j] = dydtn[j];

	Vector3d ST_F;
	ST_F = this->ST.F + this->LD.F + this->ST.get_mg();

	Vector3d ST_T;
	ST_T = this->ST.T + this->LD.T;

	double dydts[13];
	this->ST.get_dydt_(ST_F, ST_T, dvdt, dwdt, dydts);

	for (int j = 0; j < 13; j++)
		dydt[j + 13] = dydts[j];

	return;
}

//void BS_BallScrew::Shaft_Lock(void) {
//
//	this->ST.set_const(true, true);
//
//	return;
//};

// 外部荷重を代入
void BS_BallScrew::set_load(double *F, double *T) {
	for (int i = 0; i < 3; i++) {
		this->LD.F[i] = F[i];
		this->LD.T[i] = T[i];
	}
	return;
};

// 現在の状態を外部変数に入力．
void BS_BallScrew::save(BS_Out&OUT) {

	for (int i = 0; i < this->nCC; i++)
		for (int j = 0; j < this->CC[i].nBL; j++)
			this->CC[i].BL[j].save(OUT.CC[i].BL[j].x, OUT.CC[i].BL[j].v, OUT.CC[i].BL[j].q, OUT.CC[i].BL[j].w, OUT.CC[i].BL[j].ax, OUT.CC[i].BL[j].F, OUT.CC[i].BL[j].T);

	this->NT->save(OUT.NT.x, OUT.NT.v, OUT.NT.q, OUT.NT.w, OUT.NT.ax, OUT.NT.F, OUT.NT.T);

	for (int i = 0; i < this->NT->nCY; i++)
		this->NT->CY[i].save(OUT.NT_CY[i]);

	this->ST.save(OUT.ST.x, OUT.ST.v, OUT.ST.q, OUT.ST.w, OUT.ST.ax, OUT.ST.F, OUT.ST.T);

	for (int i = 0; i < this->ST.nCY; i++)
		this->ST.CY[i].save(OUT.ST_CY[i]);

	for (int i = 0; i < this->nP; i++)
		this->BNP[i].save(OUT.BNP[i]);

	for (int i = 0; i < this->nP; i++)
		this->BSP[i].save(OUT.BSP[i]);

	return;
}

BS_BallScrew::BS_BallScrew() {
	this->CC = NULL;
	this->BBP = NULL;
	this->BNP = NULL;
	this->BSP = NULL;
	return;
}

BS_BallScrew::~BS_BallScrew() {

	if (this->CC != NULL)
		delete[] this->CC;
	if (this->BBP != NULL)
		delete[] this->BBP;
	if (this->BNP != NULL)
		delete[] this->BNP;
	if (this->BSP != NULL)
		delete[] this->BSP;
	return;
}













//// 遠心力とジャイロモーメント．とりあえずベタ書きで．
//Vector3d x = this->BL[ib].x;
//Vector3d v = this->BL[ib].v;
//Vector3d w = this->BL[ib].w;
//double rv = x[1] * v[2] - x[2] * v[1];
//Vector3d x_ = Vector3d(0.0, x[1], x[2]);
//Vector3d e = x_.normalized();
//double r = x_.norm();
//double Fc = this->BL[ib].m * rv * rv / r / r / r;
//Vector3d W = Vector3d(rv / r / r, 0.0, 0.0);
//Vector3d L = this->BL[ib].I.cwiseProduct(w);
//Vector3d Tc = W.cross(L);

//// stepp（与圧計算）用のインタフェイス．
//void BS_BallScrew::get_yp(double*yp) {
//
//	for (int i = 0; i < 2; i++)
//		yp[i] = this->NT->CY[i].x.x();
//
//	for (int i = 0; i < this->nBL; i++) {
//		int i2 = i * 2;
//		Vector3d eta = this->BSP[i].get_eta();
//		yp[i2 + 2] = eta[1];
//		yp[i2 + 3] = eta[2];
//	}
//	return;
//}
//
//// stepp（与圧計算）用のインタフェイス．
//void BS_BallScrew::set_yp(const double*yp) {
//
//	for (int i = 0; i < 2; i++)
//		this->NT->dx[i] = Vector3d(yp[i], 0.0, 0.0);
//	this->NT->set_dx();
//
//	for (int i = 0; i < this->nBL; i++) {
//		int i2 = i * 2;
//		Vector2d eta = Vector2d(yp[i2 + 2], yp[i2 + 3]);
//		this->BSP[i].set_etap(eta);
//	}
//	return;
//}
//
//// 荷重を配列にして返すメソッド．
//void BS_BallScrew::get_Fp(double*Fp) {
//
//	this->NT->CY[0].F = Vector3d(-this->F_pre, 0.0, 0.0);
//	this->NT->CY[1].F = Vector3d(this->F_pre, 0.0, 0.0);
//
//	for (int i = 0; i < this->nBL; i++) {
//
//		Vector2d Fbn, Fbs;
//		Vector3d Fnb, Tnb, Fsb, Tsb;
//
//		this->BNP[i].get_F0(Fbn, Fnb, Tnb);
//		this->BSP[i].get_F0(Fbs, Fsb, Tsb);
//
//		Vector2d Fb = Fbn + Fbs;
//
//		int i2 = i * 2;
//		Fp[i2 + 2] = Fb[0];
//		Fp[i2 + 3] = Fb[1];
//
//		this->BNP[i].CY->F += Fnb;
//	}
//	for (int i = 0; i < 2; i++)
//		Fp[i] = this->NT->CY[i].F.x();
//
//	return;
//}

//
//for (size_t i = 0; i < this->nBL - 1; i++) {
//	Vector3d hoge = this->BNP[i].get_e();
//	Vector3d fuga = this->BL[i+1].x - this->BL[i].x;
//	cout << hoge.dot(fuga) / fuga.norm() << endl;
//}
//
	//for (int i = 0; i < this->nBL; i++) {
	//	Vector3d e = this->BNP[i].get_e();
	//	this->BNP[i].set_e(e);
	//	this->BSP[i].set_e(e);
	//}


//// step1（剛性計算）用のインタフェイス．
//void BS_BallScrew::get_y1(double*y1) {
//
//	int ib = this->ib;
//
//	Vector3d eta = this->BNP[ib].get_eta();
//	y1[0] = eta[1] / Rigid::l;
//	y1[1] = eta[2] / Rigid::l;
//
//	return;
//}
//
//void BS_BallScrew::set_y1(const double*y1) {
//
//	int ib = this->ib;
//
//	Vector2d eta = Vector2d(
//		y1[0] * Rigid::l,
//		y1[1] * Rigid::l
//	);
//	this->BNP[ib].set_eta0(eta);
//
//	return;
//}


//// step2（剛性・摩擦計算）用のインタフェイス．
//void BS_BallScrew::set_y2(const double*y2) {
//
//	int ib = this->ib;
//
//	Vector2d eta = Vector2d(
//		y2[0] * Rigid::l,
//		y2[1] * Rigid::l
//	);
//	this->BNP[ib].set_eta0(eta);
//
//	Vector3d v = Vector3d(y2[2], y2[3], y2[4]) / Rigid::t * Rigid::l;
//
//	this->BL[ib].v = v;
//
//	Vector3d w = Vector3d(y2[5], y2[6], y2[7]) / Rigid::t;
//	this->BL[ib].w = w;
//
//	return;
//}
//
//// 荷重を配列にして返すメソッド．
//void BS_BallScrew::get_F2(double*F2) {
//
//	int ib = this->ib;
//
//	Vector3d Fn, Tn, Fs, Ts, T_n, T_s;
//	this->BNP[ib].get_FT(Fn, Tn, T_n);
//	this->BSP[ib].get_FT(Fs, Ts, T_s);
//
//	// 遠心力とジャイロモーメント．とりあえずベタ書きで．
//	Vector3d x = this->BL[ib].x;
//	Vector3d v = this->BL[ib].v;
//	Vector3d w = this->BL[ib].w;
//	double rv = x[1] * v[2] - x[2] * v[1];
//	Vector3d x_ = Vector3d(0.0, x[1], x[2]);
//	Vector3d e = x_.normalized();
//	double r = x_.norm();
//	double Fc = this->BL[ib].m * rv * rv / r / r / r;
//	Vector3d W = Vector3d(rv / r / r, 0.0, 0.0);
//	Vector3d L = this->BL[ib].I.cwiseProduct(w);
//	Vector3d Tc = W.cross(L);
//
//	Vector3d Fb = Fn + Fs + e * Fc;
//	Vector3d Tb = (Tn + Ts + Tc) / Rigid::l;
//
//	F2[0] = Fb[0];
//	F2[1] = Fb[1];
//	F2[2] = Fb[2];
//	F2[3] = Tb[0];
//	F2[4] = Tb[1];
//	F2[5] = Tb[2];
//
//	F2[6] = 0.0;
//	F2[7] = 0.0;
//
//	return;
//}

//// step2（剛性・摩擦計算）用のインタフェイス．
//void BS_BallScrew::get_y2(double*y2) {
//
//	int ib = this->ib;
//
//	Vector3d eta = this->BNP[ib].get_eta();
//	y2[0] = eta[1] / Rigid::l;
//	y2[1] = eta[2] / Rigid::l;
//
//	//y2[0] = 0.011e-3;
//	//y2[1] = 0.020e-3;
//
//	Vector3d v = this->BL[ib].v * Rigid::t / Rigid::l;
//	y2[2] = v[0];
//	y2[3] = v[1];
//	y2[4] = v[2];
//
//	Vector3d w = this->BL[ib].w * Rigid::t;
//	y2[5] = w[0];
//	y2[6] = w[1];
//	y2[7] = w[2];
//
//	return;
//}
