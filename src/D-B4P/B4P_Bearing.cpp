/*******************************************************************************
!								"B4P_Bearing.cpp"
!													2018/12/18	[Core-T]	楠崎
!
!	4点接触玉軸受そのもののオブジェクト．
!
!	現在の状態での相互荷重やグローバル変数ライクな値の保持など，
!	結構な責務を負っている大変な人．（神クラスにならないように注意．）
!	特筆すべき点としてB4P_Inオブジェクトをバラして値を入れる責務はこいつ1人で負っている．
!
!*******************************************************************************/

#include "B4P_Bearing.h"

// 計算の入力条件を読み込む．外部が必要なデータは渡してやる．
void B4P_Bearing::init(const B4P_In&FI, double*x, double*ax) {

	// 各部材のジオメトリと材料の設定．
	double Cage_rmg0[3];
	this->init_parts(FI, Cage_rmg0);

	// 初期位置，初期速度を代入する．
	if (true)
		this->init_position(FI.balldia, FI.ballnum, FI.ballpcd, FI.cos_alp0, FI.omegair, FI.omegaor, Cage_rmg0);

	// 入力ファイルから初期位置を決定する．(未実装)
	else
		double hoge;

	//	変位・速度の自由度と力・モーメントの自由度の分だけ確保．
	this->nX = 13 * FI.ballnum + 26;

	// 内輪の拘束条件を設定．
	this->IR.set_const(FI.bound.v_const, FI.bound.w_const);
	// 内輪の初期値を設定
	this->IR.x = Vector3d(x);
	this->IR.set_ax(Vector3d(ax));
	double w = FI.omegair;
	this->IR.w = w * this->IR.get_ax();

	return;
}

void B4P_Bearing::allocate_Pair(int Z) {

	this->BOP = new B4P_BallOuterRingPair[Z];
	this->BIP = new B4P_BallInnerRingPair[Z];
	this->BCP = new B4P_BallCagePair[Z];

	return;
}

// B4P_Inから必要な値をb4pに入れるメソッド．ジオメトリと材料．
void B4P_Bearing::init_parts(const B4P_In&FI, double*Cage_rmg0) {

	// ボール個数だけ取り出し，配列の動的確保．
	this->Z = FI.ballnum;
	this->BL = new Ball[Z];

	// 各部材のジオメトリの設定．
	this->OR.init(FI.OR, FI.ballpcd);
	this->IR.init(FI.IR, FI.ballpcd);
	switch (FI.cage_type) {
	case B4P_In::snap_cage:
		this->CG = new bal_SnapCage;
		for (int i = 0; i < 3; i++)
			Cage_rmg0[i] = FI.Cage.rmg0[i];
		break;
	default:
		this->CG = new bal_SnapCage;
		break;
	}
	this->CG->init(FI);

	bool v_const[3], w_const[3];
	for (int i = 0; i < 3; i++) {
		v_const[i] = false;
		w_const[i] = false;
	}

	for (int i = 0; i < this->Z; i++)
		this->BL[i].init(FI.BL[i].dia, FI.BL[i].E, FI.BL[i].por, FI.BL[i].den, FI.BL[i].rms, v_const, w_const);

	this->allocate_Pair(this->Z);

	// ball_X_pairの相手を紐づける．
	for (int i = 0; i < this->Z; i++) {
		this->BIP[i].link(&this->BL[i], &this->IR);
		this->BOP[i].link(&this->BL[i], &this->OR);
		this->BCP[i].link(&this->BL[i], this->CG, i);
	}

	// ball_X_pairの初期化を行う．
	for (int i = 0; i < this->Z; i++) {
		this->BIP[i].init(FI);
		this->BOP[i].init(FI);
		this->BCP[i].init(FI);
	}

	// 他に必要な情報をb4pに保持する．
	this->pcd = FI.ballpcd;

	Vector3d Fin = Vector3d(FI.LoadIn[0], FI.LoadIn[1], FI.LoadIn[2]);
	Vector3d Nin = Vector3d(FI.LoadIn[3], FI.LoadIn[4], FI.LoadIn[5]);

	// 代入．
	this->F_load = Fin;
	this->T_load = Nin;

	// 各部材に荷重ベクトルを設定．
	Rigid::g = Vector3d(FI.rigid.g);

	return;
}


// B4P_Inから初期変位・速度を算出するメソッド．計算は中立位置，理論公転数，自転数とする．
void B4P_Bearing::init_position(
	double balldia,			// 玉基本径 
	int	   ballnum,			// 玉個数
	double ballpcd,			// 玉pcd
	double cos_alp0,		// 初期接触角の余弦[-] 
	double omegair,			// 内輪回転速度[rad/s]
	double omegaor,			// 外輪回転速度[rad/s]
	const double*rmg0		// 保持器重心を基準に保持器幾何中心へ向かうベクトル
) {
	// 初期クォータニオンは全て [1,0,0,0] とする．
	Quaterniond q0 = Quaterniond::Identity();

	// 位置の算出．玉はPCD上に並ぶものとする．0番目の玉は12時（+z軸）にあるものとする．
	double x; ArrayXd th(ballnum + 1), y(ballnum + 1), z(ballnum + 1);
	th = ArrayXd::LinSpaced(ballnum + 1, 0, 2 * Numeric::pi);
	x = 0;
	y = 0.5 * ballpcd *-th.sin();
	z = 0.5 * ballpcd * th.cos();

	// 初期方位角を設定する．
	this->azimuth0 = new double[ballnum];
	for (int i = 0; i < ballnum; i++)
		this->azimuth0[i] = th[i];

	// 速度の算出．保持器の理論公転数とする．テクニカルレポートp.248．玉は接触点間距離を直径にもつ小球と仮定する．
	double nc, Dw, vx; ArrayXd vy(ballnum + 1), vz(ballnum + 1);
	Dw = balldia * 0.5 * cos_alp0;	//「接触点間距離を直径にもつ小球」の直径

	nc = (1 - Dw / ballpcd) * 0.5 * omegair
		+ (1 + Dw / ballpcd) * 0.5 * omegaor; // 理論公転数 [rad/s]
	vx = 0;
	vy = 0.5 * ballpcd * nc *-th.cos();
	vz = 0.5 * ballpcd * nc *-th.sin();

	// 角速度の算出．理論自転数とする．テクニカルレポートp.248．玉は接触点間距離をもつ小球と仮定する．
	double na, wx, wy, wz;
	na = (ballpcd / Dw - Dw / ballpcd) * 0.5 * (omegaor - omegair); // 理論自転数 [rad/s]
	wx = na;
	wy = 0;
	wz = 0;

	Vector3d x0, v0, w0;
	// 各転動体に値を代入．
	for (int i = 0; i < this->Z; i++) {
		x0(0) = x;		x0(1) = y[i];	x0(2) = z[i];
		v0(0) = vx;		v0(1) = vy[i];	v0(2) = vz[i];
		w0(0) = wx;		w0(1) = wy;		w0(2) = wz;

		this->BL[i].set_param(x0, v0, q0, w0);
	}

	// 内外輪・保持器の初期変位・速度・姿勢は0．
	x0 *= 0;
	v0 *= 0;
	w0 *= 0;

	w0(0) = omegair;	// 角速度はx成分だけ入力値．他は0．
	this->IR.set_param(x0, v0, q0, w0);

	w0(0) = omegaor;	// 角速度はx成分だけ入力値．他は0．
	this->OR.set_param(x0, v0, q0, w0);

	x0 = - Vector3d(rmg0);	// 初期位置の重心は偏った場所に存在する．幾何中心と重心の差は入力値であるためそれを代入．
	w0(0) = nc;				// 角速度はx成分だけ理論公転数．他は0．
	this->CG->set_param(x0, v0, q0, w0);

	return;
}


// 外部変数から軸受各部材の変位・運動量を書き換える．
void B4P_Bearing::set_y(const double*y) {
	Vector3d    xi(y[0], y[1], y[2]);
	Vector3d    vi(y[3], y[4], y[5]);
	Quaterniond qi(y[6], y[7], y[8], y[9]);
	Vector3d    wi(y[10], y[11], y[12]);
	this->IR.set_y(xi, vi, qi, wi);
	//cout  << endl;
	//cout << "*qi=**************************" << endl;
	//std::cout << IR.q.x() << std::endl;
	//std::cout << IR.q.y() << std::endl;
	//std::cout << IR.q.z() << std::endl;
	//std::cout << IR.q.w() << std::endl;
	//cout << "wi=**************************" << endl;
	//
	//std::cout << IR.w << std::endl;
	//cout << "ax=**************************" << endl;
	//std::cout << IR.get_ax() << std::endl;

	Vector3d    xc(y[13], y[14], y[15]);
	Vector3d    vc(y[16], y[17], y[18]);
	Quaterniond qc(y[19], y[20], y[21], y[22]);
	Vector3d    wc(y[23], y[24], y[25]);
	this->CG->set_y(xc, vc, qc, wc);

	Vector3d xb, vb, wb;
	Quaterniond qb;
	for (int i = 0; i < this->Z; i++) {
		int j = 13 * i + 26;
		xb = Vector3d(y[j + 0], y[j + 1], y[j + 2]);
		vb = Vector3d(y[j + 3], y[j + 4], y[j + 5]);
		qb = Quaterniond(y[j + 6], y[j + 7], y[j + 8], y[j + 9]);
		wb = Vector3d(y[j + 10], y[j + 11], y[j + 12]);
		this->BL[i].set_y(xb, vb, qb, wb);
	}
	return;
}

// 軸受各部材の変位・運動量を外部変数に出力する．
void B4P_Bearing::get_y(double*y) {

	Vector3d x, p, w; Quaterniond q;

	for (int i = 0; i < this->Z; i++) {
		this->BL[i].get_y(x, p, q, w);
		y[13 * i + 26] = x.x(); y[13 * i + 27] = x.y(); y[13 * i + 28] = x.z();
		y[13 * i + 29] = p.x(); y[13 * i + 30] = p.y(); y[13 * i + 31] = p.z();
		y[13 * i + 32] = q.w(); y[13 * i + 33] = q.x(); y[13 * i + 34] = q.y(); y[13 * i + 35] = q.z();
		y[13 * i + 36] = w.x(); y[13 * i + 37] = w.y(); y[13 * i + 38] = w.z();
	}
	this->IR.get_y(x, p, q, w);
	y[0] = x.x(); y[1] = x.y(); y[2] = x.z();
	y[3] = p.x(); y[4] = p.y(); y[5] = p.z();
	y[6] = q.w(); y[7] = q.x(); y[8] = q.y(); y[9] = q.z();
	y[10] = w.x(); y[11] = w.y(); y[12] = w.z();

	this->CG->get_y(x, p, q, w);
	y[13] = x.x(); y[14] = x.y(); y[15] = x.z();
	y[16] = p.x(); y[17] = p.y(); y[18] = p.z();
	y[19] = q.w(); y[20] = q.x(); y[21] = q.y(); y[22] = q.z();
	y[23] = w.x(); y[24] = w.y(); y[25] = w.z();

	return;
}

// 現在の状態での各部材にかかる荷重を計算し，微分値を返すメソッド．
void B4P_Bearing::get_dydt(double*dydt) {

	Vector3d sumFo, sumTo, sumFi, sumTi, sumFc, sumTc;
	sumFo = sumTo = sumFi = sumTi = sumFc = sumTc = Vector3d::Zero();

	for (int i = 0; i < this->Z; i++) {

		// 各ボールの荷重と摩擦を計算．
		Vector3d Fbi, Tbi, Fib, Tib; // Fbi : ボール(b)が内輪(i)から受ける力．
		this->BIP[i].calc_force(Fbi, Tbi, Fib, Tib);

		Vector3d Fbo, Tbo, Fob, Tob; // Fbo : ボール(b)が外輪(o)から受ける力．
		this->BOP[i].calc_force(Fbo, Tbo, Fob, Tob);

		Vector3d Fbc, Tbc, Fcb, Tcb; // Fbc : ボール(b)が保持器(c)から受ける力．
		this->BCP[i].calc_force(Fbc, Tbc, Fcb, Tcb);

		// ボールにかかる荷重の総和を算出．
		Vector3d sumFb = Fbi + Fbo + Fbc;
		Vector3d sumTb = Tbi + Tbo + Tbc;

		// 重力加速度を加算．
		sumFb += this->BL[i].get_mg();

		double dydt_[13];
		this->BL[i].get_dydt(sumFb, sumTb, dydt_);

		// ループごとに戻り値に代入．
		for (int j = 0; j < 13; j++)
			dydt[13 * i + 26 + j] = dydt_[j];

		// 各部材にかかる荷重を積算．
		sumFo += Fob;	sumTo += Tob;
		sumFi += Fib;	sumTi += Tib;
		sumFc += Fcb;	sumTc += Tcb;

		// 書き込み用の変数に保存．
		this->BL[i].set_FT(sumFb, sumTb);
	}

	// 重力加速度を加算．
	sumFo += this->OR.get_mg();
	sumFi += this->IR.get_mg();
	sumFc += this->CG->get_mg();

	// 内輪だけ入力荷重を加算．（いずれ外輪・保持器にも実装予定？）
	sumFi += this->F_load;
	sumTi += this->T_load;

	double dyidt_[13];
	this->IR.get_dydt_(sumFi, sumTi, dyidt_);

	for (int j = 0; j < 13; j++)
		dydt[0 + j] = dyidt_[j];

	double dycdt_[13];
	this->CG->get_dydt(sumFc, sumTc, dycdt_);

	for (int j = 0; j < 13; j++)
		dydt[13 + j] = dycdt_[j];

	// 書き込み用の変数に保存．
	this->OR.set_FT(sumFo, sumTo);
	this->IR.set_FT(sumFi, sumTi);
	this->CG->set_FT(sumFc, sumTc);

	return;
}

// 外部に書き込むべき変数を Eigen の配列にして外に渡すメソッド．
void B4P_Bearing::save(B4P_Out&OUT) {

	for (int i = 0; i < this->Z; i++)
		this->BL[i].save(OUT.BL[i].x, OUT.BL[i].v, OUT.BL[i].q, OUT.BL[i].w, OUT.BL[i].ax, OUT.BL[i].F, OUT.BL[i].T);

	this->OR.save(OUT.OR.x, OUT.OR.v, OUT.OR.q, OUT.OR.w, OUT.OR.ax, OUT.OR.F, OUT.OR.T);
	this->IR.save(OUT.IR.x, OUT.IR.v, OUT.IR.q, OUT.IR.w, OUT.IR.ax, OUT.IR.F, OUT.IR.T);
	this->CG->save(OUT.CG.x, OUT.CG.v, OUT.CG.q, OUT.CG.w, OUT.CG.ax, OUT.CG.F, OUT.CG.T);
	
	for (int i = 0; i < this->Z; i++) {
		this->BOP[i].write(OUT.BOP[i]);
		this->BIP[i].write(OUT.BIP[i]);
		this->BCP[i].save(OUT.BCP[i]);
	}
	return;
}

// 静解析(1) 剛性計算：軸受各部材の変位・運動量を外部変数に出力
void B4P_Bearing::get_Xstf(double*Xstf) {
	this->IR.get_Xstf(Xstf);

	for (int i = 0; i < this->Z; i++) {
		Vector3d x, v, w;
		Quaterniond q;
		this->BL[i].get_param(x, v, q, w);
		Vector3d thXZ = this->OR.ine_to_XZcoord(x);
		Xstf[i * 2 + 5] = thXZ[1] / Rigid::l;
		Xstf[i * 2 + 6] = thXZ[2] / Rigid::l;
	}
	return;
}


// 静解析(1) 剛性計算：calculatorクラスからの入力 X_stf を軸受に代入
void B4P_Bearing::set_Xstf(double* Xstf) {

	this->IR.set_Xstf(Xstf);

	// 玉は η-ξ 座標系から慣性座標系に変換
	for (int i = 0; i < this->Z; i++) {
		int i2 = i * 2;
		Vector3d eta(this->azimuth0[i], Xstf[i2 + 5] * Rigid::l, Xstf[i2 + 6] * Rigid::l);
		this->BOP[i].set_eta_stf(eta);
	}

	return;
}
// 静解析(1) 剛性計算：軸受内部荷重を計算し，calculatorクラスに F(X) の形で出力
void B4P_Bearing::get_Fstf(double* Fstf) {

	Vector3d sumFi, sumTi;
	sumFi << 0, 0, 0;	sumTi << 0, 0, 0;

	// 入力荷重と重力加速度を加算．
	sumFi += this->F_load + this->IR.get_mg();
	sumTi += this->T_load;

	for (int i = 0; i < this->Z; i++) {

		// 各ボールの荷重を計算．
		Vector3d Fbi, Fbo;				// Fbi : ボール(b)が内輪(i)から受ける力（玉軌道座標系）
		Vector3d Fib, Tib, Fob, Tob;	// Fib : 内輪(i)がボール(b)から受ける力（慣性座標系）
		this->BIP[i].get_Fstf(Fbi, Fib, Tib);
		this->BOP[i].get_Fstf(Fbo, Fob, Tob);


		// 外輪・内輪それぞれの溝直角断面で見たときの荷重を計算．
		Vector2d Fbo_ = this->OR.to_cross_sction(Fbo, this->azimuth0[i]);
		Vector2d Fbi_ = this->IR.to_cross_sction(Fbi, this->azimuth0[i]);
		Vector2d Fb = Fbo_ + Fbi_;

		// 弱いバネ（斥力）．x^-1 の形にすることによって，外側に収束しやすくする．
		double x = Vector2d(this->BL[i].x[1], this->BL[i].x[2]).norm();
		double k = 1e-20;
		Fstf[2 * i + 5] = Fb[0] - k / x;
		Fstf[2 * i + 6] = Fb[1];

		// 接触荷重計算で出てきた内輪側の反作用を積算していく．
		sumFi += Fib;
		sumTi += Tib;
	}



	// 全てのボールの荷重の総和を戻り値として与える．
	Fstf[0] = sumFi[0];
	Fstf[1] = sumFi[1];
	Fstf[2] = sumFi[2];

	// トルクを内輪座標系に変換し，第1,2成分を戻り値にする．（第0成分は回す方向なので考慮しない）
	Vector3d Ti = this->IR.to_myvector(sumTi);
	Fstf[3] = Ti[1];
	Fstf[4] = Ti[2];
	return;
}

B4P_Bearing::B4P_Bearing(void) {
	this->BL = NULL;
	this->BOP = NULL;
	this->BIP = NULL;
	this->BCP = NULL;
	this->azimuth0 = NULL;
	return;
}

B4P_Bearing::~B4P_Bearing(void) {
	if (this->BL != NULL)
		delete[] this->BL;
	if (this->BOP != NULL)
		delete[] this->BOP;
	if (this->BIP != NULL)
		delete[] this->BIP;
	if (this->BCP != NULL)
		delete[] this->BCP;
	if (this->azimuth0 != NULL)
		delete[] this->azimuth0;
	return;
}































//this->nX_stf =  2 * FI.ballnum + 5;
//this->nX_frc =  4 * FI.ballnum;
//this->nX_stf_frc =  6 * FI.ballnum + 5;


//// YZ方向外部荷重の単位ベクトルを求める（1:Y方向，2:Z方向）
//Vector3d B4P_Bearing::get_Fyzload_dir() {
//	// YZ平面上の力の大きさ
//	double F_str = sqrt(this->F_load[1]*this->F_load[1]+this->F_load[2]*this->F_load[2]);
//	if (F_str == 0) return Vector3d(0, 0, 0);
//	// 単位ベクトルの導出
//	return Vector3d(0, this->F_load[1]/F_str, this->F_load[2]/F_str);
//}
//
//// X方向外部荷重の向きを求める（-1:負方向，1:正方向，0:X荷重無）
//int B4P_Bearing::get_Fxload_dir() {
//
//	if (this->F_load[0] > 0) return 1;
//	else if (this->F_load[0] < 0) return -1;
//	else return 0;
//}
//
//// 近似解3.YZ方向外部荷重モーメントの方向ベクトルを求める（1:My方向，2:Mz方向）
//Vector3d B4P_Bearing::get_nyzload_dir() {
//	// YZ平面上の力の大きさ
//	double N_str = sqrt(this->T_load[1]*this->T_load[1]+this->T_load[2]*this->T_load[2]);
//	if (N_str == 0) return Vector3d(0, 0, 0);
//	// 単位ベクトルの導出
//	return Vector3d(0, this->T_load[1]/N_str, this->T_load[2]/N_str);
//}
//
//
//
//// 近似解3. 軸方向の姿勢を入力
//void B4P_Bearing::set_q_ir(Vector3d ax) {
//	this->IR.set_ax(ax);
//	return;
//}

//
//
//// 現在の状態での内輪・玉１つにかかる荷重を算出するメソッド．（全て慣性座標系．F[0~2]：内輪荷重，F[3~4]：内輪トルク，F[5~6]：玉荷重．）
//void B4P_Bearing::get_F_stf_ball(double*F, int i) {
//
//	// 各ボールの荷重を計算．
//	Vector3d Fbi, Fib, Tib; // Fbi : ボール(b)が内輪(i)から受ける力．
//	this->BIP[i].calc_force_stf(Fbi, Fib, Tib);
//	Vector3d Fbo, Fob, Tob; // Fbo : ボール(b)が外輪(o)から受ける力．
//	this->BOP[i].calc_force_stf(Fbo, Fob, Tob);
//
//	// 重力加速度を加算．いちおう外輪基準ということで，外輪荷重に加算している．
//	// さらに遠心力を加算．
//	Fbo += this->BL[i].get_mg() + this->BOP[i].get_centf();
//
//	// 外輪・内輪それぞれの溝直角断面で見たときの荷重を計算．
//	Vector2d Fbo_ = this->OR.to_cross_sction(Fbo, this->azimuth[i]);
//	Vector2d Fbi_ = this->IR.to_cross_sction(Fbi, this->azimuth[i]);
//	Vector2d Fb   = Fbo_ + Fbi_;
//
//
//
//
//	// ボールにかかる荷重の総和を算出．
//	F[5] = Fb[0];
//	F[6] = Fb[1];
//
//	// 全てのボールの荷重の総和を戻り値として与える．
//	F[0] = Fib[0];
//	F[1] = Fib[1];
//	F[2] = Fib[2];
//
//	// トルクを内輪座標系に変換し，第1,2成分を戻り値にする．（第0成分は回す方向なので考慮しない）
//	Vector3d Ti = this->IR.to_myvector(Tib);
//	F[3] = Ti[1];
//	F[4] = Ti[2];
//}
//
//
//// 外部変数から軸受各部材の変位・運動量を書き換える．（剛性収束計算：STEP1用．）
//void B4P_Bearing::set_param_stf(const double*y_stf) {
//	//内輪の各変数を書き換える
//	Vector3d    xi(y_stf[0], y_stf[1], y_stf[2]);
//	Vector3d    v0(0.0, 0.0, 0.0);
//	Quaterniond qi(1.0, 0.0, y_stf[3], y_stf[4]);
//	//Vector3d    w0( 1000.0,  0.0,  0.0);
//	Vector3d    w0(0.0, 0.0, 0.0);
//
//	this->IR.set_param(xi, v0, qi, w0);
//
//	for (int i = 0; i < this->Z; i++) {
//		Vector3d xb = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5], y_stf[2*i+6]));
//		this->BL[i].set_param_stf(xb);
//	}
//}
//
//
//// 変位が y_stf であるときの Jacobian を算出するメソッド．
//void B4P_Bearing::get_Jacobian_stf(double*y_stf, double dx, double dq, double*Jacobian) {
//
//	int N = this->nX_stf;
//
//	for (int i = 0; i < N * N; i++)
//		Jacobian[i] = 0.0;
//
//	double F_Xp[7], F_Xm[7], F_Zp[7], F_Zm[7];
//	for (int i = 0; i < this->Z; i++) {
//		Vector3d x = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5], y_stf[2*i+6]));
//		this->BL[i].set_param_stf(x);
//	}
//
//
//	for (int i = 0; i < this->Z; i++) {
//		// X方向に微小変位+dxを追加
//		Vector3d Xp = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5]+dx, y_stf[2*i+6]));
//		this->BL[i].set_param_stf(Xp);
//		this->get_F_stf_ball(F_Xp, i);
//
//		// X方向に微小変位-dxを追加
//		Vector3d Xm = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5]-dx, y_stf[2*i+6]));
//		this->BL[i].set_param_stf(Xm);
//		this->get_F_stf_ball(F_Xm, i);
//
//		// Xについて偏微分をとる．(F[0~2]：内輪荷重，F[3~4]：内輪トルク，F[5~6]：玉荷重，ということを考慮に入れて計算する）
//		//for (int j = 0; j < 5; j++){
//		int j;
//		j = 0;	Jacobian[N * j + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		j = 1;	Jacobian[N * j + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		j = 2;	Jacobian[N * j + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		j = 3;	Jacobian[N * j + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		j = 4;	Jacobian[N * j + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		j = 5;	Jacobian[N * (2 * i + 5) + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		j = 6;	Jacobian[N * (2 * i + 6) + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		//}
//
//		// Z方向に微小変位+dxを追加
//		Vector3d Zp = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5], y_stf[2*i+6]+dx));
//		this->BL[i].set_param_stf(Zp);
//		this->get_F_stf_ball(F_Zp, i);
//
//		// Z方向に微小変位-dxを追加
//		Vector3d Zm = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5], y_stf[2*i+6]-dx));
//		this->BL[i].set_param_stf(Zm);
//		this->get_F_stf_ball(F_Zm, i);
//
//		// Zについて偏微分をとる．
//		//int j;
//		j = 0;	Jacobian[N * j + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//		j = 1;	Jacobian[N * j + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//		j = 2;	Jacobian[N * j + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//		j = 3;	Jacobian[N * j + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//		j = 4;	Jacobian[N * j + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//		j = 5;	Jacobian[N * (2 * i + 5) + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//		j = 6;	Jacobian[N * (2 * i + 6) + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//
//
//		// ボールを元の位置に戻す．
//		Vector3d x = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5], y_stf[2*i+6]));
//		this->BL[i].set_param_stf(x);
//	}
//
//	Vector3d    xi(y_stf[0], y_stf[1], y_stf[2]);
//	Vector3d    v0(0.0, 0.0, 0.0);
//	Quaterniond qi(1.0, 0.0, y_stf[3], y_stf[4]);
//	Vector3d    w0(0.0, 0.0, 0.0);
//	double*F_p = new double[N];
//	double*F_m = new double[N];
//
//	for (int i = 0; i < 3; i++) {
//		xi[i] = y_stf[i] + dx;
//		this->IR.set_param(xi, v0, qi, w0);
//		this->get_F_stf(F_p);
//
//		xi[i] = y_stf[i] - dx;
//		this->IR.set_param(xi, v0, qi, w0);
//		this->get_F_stf(F_m);
//
//		xi[i] = y_stf[i];
//		this->IR.set_param(xi, v0, qi, w0);
//
//		for (int j = 0; j < N; j++) {
//			Jacobian[N * i + j] = (F_p[j] - F_m[j]) / (2 * dx);
//		}
//	}
//
//	// クォータニオンはy,zそれぞれ計算．
//	qi.y() = y_stf[3] + dx;
//	this->IR.set_param(xi, v0, qi, w0);
//	this->get_F_stf(F_p);
//
//	qi.y() = y_stf[3] - dx;
//	this->IR.set_param(xi, v0, qi, w0);
//	this->get_F_stf(F_m);
//
//	for (int j = 0; j < N; j++) {
//		Jacobian[N * 3 + j] = (F_p[j] - F_m[j]) / (2 * dx);
//		//cout << N * 3 + j << ":" <<
//		//	Jacobian[N * 3 + j] << endl;
//	}
//	qi.y() = y_stf[3];
//	qi.z() = y_stf[4] + dx;
//	this->IR.set_param(xi, v0, qi, w0);
//	this->get_F_stf(F_p);
//
//	qi.z() = y_stf[4] - dx;
//	this->IR.set_param(xi, v0, qi, w0);
//	this->get_F_stf(F_m);
//
//	for (int j = 0; j < N; j++) {
//		Jacobian[N * 4 + j] = (F_p[j] - F_m[j]) / (2 * dx);
//		//cout << N * 4 + j << ":" <<
//		//	Jacobian[N * 4 + j] << endl;
//	}
//	// 元に戻す
//	this->set_param_stf(y_stf);
//
//}
//
//



//// 外部変数から軸受各部材の変位・運動量を書き換える．（摩擦収束計算：STEP2用．）
//void B4P_Bearing::set_param_frc(const double*y_frc) {
//
//	for (int i = 0; i < this->Z; i++) {
//		int j = 4 * i;
//		Vector3d vb = this->OR.Y_to_inevelocity(y_frc[j+0], this->azimuth[i]);
//		Vector3d wb = Vector3d(y_frc[j+1], y_frc[j+2], y_frc[j+3]);
//		this->BL[i].set_param_frc(vb, wb);
//	}
//}
//
//// 保持器公転数を取得するメソッド．
//double B4P_Bearing::get_wCage(void) {
//	return -1e5;
//}
//
//// 軸受各部材の変位・運動量を外部変数に出力する．（摩擦収束計算：STEP2用．）
//void B4P_Bearing::get_param_frc(double*y_frc) {
//	Vector3d vb, wb;
//	for (int i = 0; i < this->Z; i++) {
//		int j = 4 * i;
//		this->BL[i].get_param_frc(vb, wb);
//		double v = this->OR.ine_to_Yvelocity(vb, this->azimuth[i]);
//		y_frc[j+0] = v; y_frc[j+1] = wb[0]; y_frc[j+2] = wb[1]; y_frc[j+3] = wb[2];
//	}
//}
	/*旧inputファイルからの読み込み*/
	// ボール個数だけ取り出し，配列の動的確保．
	/*
	int Z     = FI.Z;
	this->Z   = Z;
	this->BL  = new Ball[Z];
	this->BOP = new B4P_BallRingPair[Z];
	this->BIP = new B4P_BallRingPair[Z];
	this->BCP = new B4P_BallCagePair[Z];
	*/
	/*新inputファイルからの読み込み*/
