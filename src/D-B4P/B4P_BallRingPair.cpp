#include "B4P_BallRingPair.h"

// 玉・内輪（外輪）オブジェクトを所有させる
void B4P_BallRingPair::link(Ball*BL, B4P_Ring*RG) {
	this->BL = BL;
	this->RG = RG;

	return;
}

// 現在の状態での玉の接触判定および，そのときの接触状態を計算
// 参考 "01. 円環面と球体の弾性接触アルゴリズム.pptx"
bool B4P_BallRingPair::how_Contact
(							// out:	[bit]	接触していればtrue，していなければfalseを返す．
	int i,					// in:	[-]		接触溝番号．0 or 1．
	const Vector3d&bl_x,	// in:	[m]		ボール位置．（リング座標系）
	Vector3d&er,			// out:	[-]:	リング中心から見たボール位相方向ベクトル．x成分=0．(リング座標系基準)
	Vector3d&eg,			// out:	[-]:	溝中心から見たボール方向ベクトル．ボールの受ける接触荷重方向の逆向き．（リング座標系基準）
	double&dx				// out:	[-]		弾性接近量．正で接触．
) {

	// まずyz平面で見たボールの方向ベクトル er (リング座標系基準)を導出．
	double r = sqrt(Numeric::Square(bl_x[1]) + Numeric::Square(bl_x[2]));
	er = Vector3d(0.0, bl_x[1] / r, bl_x[2] / r);

	// 着目する溝中心点（リング座標系基準）を決定．
	Vector3d x_trs = this->RG->GV[i].Rr *er;
	x_trs[0] = this->RG->GV[i].Rx;

	// 溝中心点からボール中心に向かうベクトルの算出．
	Vector3d gb = bl_x - x_trs;
	double   gb_norm = gb.norm();

	// 溝中心点とボール中心が一致していたら接していないことにする．（0除算防止）
	if (gb_norm == 0)
		return false;

	eg = gb / gb_norm;
	// 幾何学からボール食い込み量の算出．
	dx = gb_norm - this->RG->GV[i].r + this->BL->r;

	// ボール中心がトーラスから出ていたら接していないことにする．
	if (gb_norm > this->RG->GV[i].r) {
		cout << "玉が内外輪の軌道から外れています．直ちに計算を停止してください．" << endl;
		return false;
	}

	// 接近量が負値の場合接触していない．dx=0のときも以降の計算で0割が発生するため，接触なしとする．
	if (dx <= 0)
		return false;

	return true;
}


// Hertzの計算から接触剛性を返すメソッド（減衰なし）
// 参考 "01. 円環面と球体の弾性接触アルゴリズム.pptx" STEP 4 ~ 7
void B4P_BallRingPair::calc_Hertz
(							// out:	[-]		なし
	int i,					// in:	[-]		接触溝番号．0 or 1．
	const Vector3d&bl_x,	// in:	[m]		ボール位置．（リング座標系）
	const Vector3d&er,		// in:	[-]:	リング中心から見たボール位相方向ベクトル．x成分=0．(リング座標系)
	const Vector3d&eg,		// in:	[-]:	溝中心から見たボール方向ベクトル．ボールの受ける接触荷重方向の逆向き．(リング座標系)
	double dx,				// in:	[-]		弾性接近量．正で接触．
	double&Rx,				// out:	[m]:	転がり方向曲率半径．
	double&Ry,				// out:	[m]:	周方向曲率半径．
	Vector3d&p,				// out:	[m]:	接触点位置（慣性座標系）．
	double&cos_alp,			// out:	[-]:	接触角のcos．
	double&sin_alp,			// out:	[-]:	接触角のsin．
	double&a,				// out:	[m]:	接触楕円長径．
	double&b,				// out:	[m]:	接触楕円短径．
	double&k				// out:	[N/x^1.5]:	非線形剛性．
) {
	// 接触角 alp の cos の算出．（ナットラジアル当たり=0度，シャフトラジアル当たり=180度）
	Vector3d ea(1.0, 0.0, 0.0);
	cos_alp = er.dot(eg);
	sin_alp = ea.dot(eg);

	// 曲率の算出．（BrewHamrock の近似式で使用）
	double rho0 = this->BL->r_inv;			// ボールの周方向曲率．
	double rho1 = this->BL->r_inv;			// ボールの径方向曲率．
	double rho2 = -cos_alp / (this->RG->GV[i].Rr + this->RG->GV[i].r * cos_alp);	// 輪の転がり方向曲率．
	double rho3 = -this->RG->GV[i].r_inv;	// 輪の周方向曲率．

	Rx = 1. / (rho0 + rho2);
	Ry = 1. / (rho1 + rho3);

	// BrewHamrock の近似式から剛性，接触楕円を求める． 
	this->HZ->calc(Rx, Ry, dx, this->E, k, a, b);

	// 比較用．以前までの仕様です．確認が終わったら以下2行は消してください．
	//double Rm = 2 * this->BL->r * this->RG->GV[i].r / (this->BL->r + this->RG->GV[i].r);	// 接触面の曲率半径
	//Vector3d p__ = bl_x + (this->BL->r + a * a * 0.5 * (1 / Rm - 1 / this->BL->r)) * eg;

	// 接触点位置pの算出．（リング座標系⇒慣性座標系．）
	double x0 = Numeric::EffectiveCenter(this->BL->r, this->RG->GV[i].r, dx);
	Vector3d p_ = bl_x + x0 * eg;
	p = this->RG->to_inecoord(p_);

	return;
}

// Hertzの計算から，ダンパによる減衰も考慮した荷重を返すメソッド．
double B4P_BallRingPair::calc_DynamicHertz
(							// out:	[N]		接触荷重（スカラー）
	int i,					// in:	[-]		接触溝番号．0 or 1．
	const Vector3d&bl_x,	// in:	[m]		ボール位置．（リング座標系）
	const Vector3d&bl_v,	// in:	[m/s]	ボール速度．（リング座標系）
	const Vector3d&er,		// in:	[-]:	リング中心から見たボール位相方向ベクトル．x成分=0．(リング座標系)
	const Vector3d&eg,		// in:	[-]:	溝中心から見たボール方向ベクトル．ボールの受ける接触荷重方向の逆向き．(リング座標系)
	double dx,				// in:	[-]		弾性接近量．正で接触．
	double&Rx,				// out:	[m]:	転がり方向曲率半径．
	double&Ry,				// out:	[m]:	周方向曲率半径．
	Vector3d&p,				// out:	[m]:	接触点位置（慣性座標系）．
	double&cos_alp,			// out:	[-]:	接触角のcos．
	double&sin_alp,			// out:	[-]:	接触角のsin．
	double&a,				// out:	[m]:	接触楕円長径．
	double&b				// out:	[m]:	接触楕円短径．
) {
	// Hertz接触剛性を計算
	double k;
	this->calc_Hertz(i, bl_x, er, eg, dx, Rx, Ry, p, cos_alp, sin_alp, a, b, k);

	// 玉速度，接触面垂直成分（リング座標系とリング座標系の内積）
	double bl_vn = bl_v.dot(eg);

	// 非線形の減衰力(スカラー)
	double F_norm_ = this->DF->calc(k, this->zeta, this->m, bl_vn, dx);
	return F_norm_;
}

// 滑り摩擦を求めるメソッド，トルクは軌道輪にかかるものと玉にかかるものを両方計算
void B4P_BallRingPair::calc_Sliding
(
	int i,				// in:	[-]		接触溝番号．0 or 1．
	const Vector3d&p,	// in:	[m]:	接触点位置（慣性座標系）．
	double a,			// in:	[m]:	接触楕円長径．
	double Pmax,		// in:	[Pa]:	最大面圧．
	double F_norm,		// in:	[N]:	ヘルツ接触荷重．
	double fratio,		// in:	[-]:	油膜接触割合．
	Vector3d&Fbs,		// out:	[N]:	ボールにかかる滑り摩擦力（慣性座標系）
	Vector3d&Tbs,		// out:	[Nm]:	ボールにかかる滑り摩擦トルク（慣性座標系）
	Vector3d&Tis		// out:	[Nm]:	内輪(外輪)にかかる滑り摩擦トルク（慣性座標系）
) {
	// 接触楕円をスライスしてそれぞれの滑り摩擦力を積算する．
	Fbs = Tbs = Tis = Vector3d::Zero();
	this->RG->calc_slice(p, a, this->GV[i].msmax, i, this->BL->r, this->GV[i].ps);

	for (int j = 0; j < this->GV[i].msmax; j++) {
		Vector3d us_, ur_;
		double f_arr = F_norm * this->GV[i].ratio_slice[j];
		this->get_us_ur(this->GV[i].ps[j], us_, ur_);
		double us_norm = us_.norm(), ur_norm = ur_.norm();
		double mu_tr = this->TR->calc(this->lb_eta, Pmax, us_norm, ur_norm);
		double mu_cl = this->CL->calc(this->mu, us_norm, ur_norm, this->clmb_s);
		double Fs_norm = mu_tr * fratio * f_arr + mu_cl * (1.0 - fratio) * f_arr;
		Vector3d Fs_ = Fs_norm * -us_ / us_norm;
		Vector3d Tbs_ = this->BL->calc_Torque(this->GV[i].ps[j], Fs_);
		Vector3d Tis_ = this->RG->calc_Torque(this->GV[i].ps[j], -Fs_);
		Fbs += Fs_;
		Tbs += Tbs_;
		Tis += Tis_;
		this->save_Slice(i, j, f_arr, Fs_, Tbs_, mu_cl, mu_tr, us_, this->GV[i].ps[j]);
	}
	return;
}

// 計算結果をメンバ変数に格納するメソッド．
void B4P_BallRingPair::save(bool c, int i, double cos_alp, double sin_alp, double dx, double a, double b, const Vector3d&p, const Vector3d&Fn, double Pmax, const Vector3d&us, const Vector3d&ur, const Vector3d&Tr, double lambda, double fratio, const Vector3d&Fs, const Vector3d&Ts, const Vector3d&Fb, const Vector3d&Tb, const Vector3d&Ti
) {
	this->GV[i].cos_alp = cos_alp;
	this->GV[i].sin_alp = sin_alp;
	this->GV[i].dx = dx;
	this->GV[i].a = a;
	this->GV[i].b = b;
	this->GV[i].p = p;
	this->GV[i].Fn = Fn;
	this->GV[i].Pmax = Pmax;
	this->GV[i].us = us;
	this->GV[i].ur = ur;
	this->GV[i].Tr = Tr;
	this->GV[i].lambda = lambda;
	this->GV[i].fratio = fratio;
	this->GV[i].Fs = Fs;
	this->GV[i].Ts = Ts;
	this->GV[i].Fb = Fb;
	this->GV[i].Tb = Tb;
	this->GV[i].Ti = Ti;

	return;
}
// スライス片の結果に0を格納
void B4P_BallRingPair::init_Sliceparam(int i) {
	for (int j = 0; j < this->GV[i].msmax; j++) {
		this->GV[i].SL[j].f_arr = 0;
		this->GV[i].SL[j].fs = Vector3d::Zero();	// 滑り摩擦力[N](慣性座標系)
		this->GV[i].SL[j].ts = Vector3d::Zero();	// 滑り摩擦トルク[Nm](慣性座標系)
		this->GV[i].SL[j].mu_cl = 0;				// クーロン摩擦係数[-]
		this->GV[i].SL[j].mu_tr = 0;				// トラクション係数[-]
		this->GV[i].SL[j].us = Vector3d::Zero();	// 滑り速度[m/s](慣性座標系)
		this->GV[i].SL[j].ps = Vector3d::Zero();	// 滑り速度[m/s](慣性座標系)			
	}
}



// スライス計算結果をメンバ変数に格納
void B4P_BallRingPair::save_Slice(
	int i,			// in: 溝番号
	int j,			// in: スライス番号
	double f_arr,	// in: 各スライスの接触荷重[N]
	const Vector3d& Fs_,	// in: 滑り摩擦荷重[N]
	const Vector3d& Ts_,	// in: 滑り摩擦トルク[Nm]
	double mu_cl,	// in: クーロン摩擦係数[-]
	double mu_tr,	// in: トラクション係数[-]
	const Vector3d& us_,	// in: 滑り速度[m/s]
	const Vector3d& ps_		// in: スライス片中心[m]
) {
	this->GV[i].SL[j].f_arr = f_arr;
	this->GV[i].SL[j].fs = Fs_;			// 滑り摩擦力[N](慣性座標系)
	this->GV[i].SL[j].ts = Ts_;			// 滑り摩擦トルク[Nm](慣性座標系)
	this->GV[i].SL[j].mu_cl = mu_cl;		// クーロン摩擦係数[-]
	this->GV[i].SL[j].mu_tr = mu_tr;		// トラクション係数[-]
	this->GV[i].SL[j].us = us_;			// 滑り速度[m/s](慣性座標系)
	this->GV[i].SL[j].ps = ps_;		// スライス片の相対位置[m](慣性座標系)
	return;
}

// 現在のボール-リングの変位から，お互いにかかる荷重を算出する．
void B4P_BallRingPair::calc_force(
	Vector3d&Fbi,	// out:	[N]:	ボールにかかる全ての力（慣性座標系）
	Vector3d&Tbi,	// out:	[Nm]:	ボールにかかる全てのトルク（慣性座標系）
	Vector3d&Fib, 	// out:	[N]:	リングにかかる全ての力（慣性座標系）
	Vector3d&Tib	// out:	[Nm]:	リングにかかる全てのトルク（慣性座標系）
) {
	Vector3d bl_x = this->RG->to_mycoord(this->BL->x);
	Vector3d bl_v = this->RG->to_myvelocity(this->BL->v);
	Fbi = Tbi = Fib = Tib = Vector3d::Zero();

	// groove と ball の接触計算．
	for (int i = 0; i < 2; i++) {

		// 必要な戻り値の確保．
		double dx, Rx, Ry, cos_alp, sin_alp, a, b, Pmax, lambda, fratio, k;
		dx = Rx = Ry = sin_alp = a = b = Pmax = lambda = fratio = k = 0;
		cos_alp = 1.0;
		Vector3d er, eg, p, Fn, ur, us, ur3, us3, Tbr, Fbs, Tbs, Tis, Fb, Tb, Ti;
		er = eg = p = Fn = ur = us = Tbr = Fbs = Tbs = Tis = Fb = Tb = Ti = Vector3d::Zero();

		// 接していているか判定．
		bool c = this->how_Contact(i, bl_x, er, eg, dx);
		if (c) {
			// 接触荷重・転動体面圧を導出（慣性座標系）
			double F_norm = this->calc_DynamicHertz(i, bl_x, bl_v, er, eg, dx, Rx, Ry, p, cos_alp, sin_alp, a, b);
			Fn = F_norm * this->RG->to_inevector(-eg);
			Pmax = 1.5 * F_norm / (Numeric::pi * a * b);

			//// 着目する溝中心点（リング座標系基準）を決定．
			//Vector3d x_trs = this->RG->GV[i].Rr *er;
			//x_trs[0] = this->RG->GV[i].Rx;
			//// 溝中心点からボール中心に向かうベクトルの算出．
			//Vector3d gb = bl_x - x_trs;
			//
			//Vector3d _bl_x = Vector3d(0, bl_x[1], bl_x[2]);
			//Vector3d _bl_v = Vector3d(0, bl_v[1], bl_v[2]);
			//Vector3d _dpdt = bl_v + this->BL->r * (_bl_v - this->RG->GV[i].Rr * _bl_v / _bl_x.norm()) / gb.norm();
			//Vector3d dpdt = this->RG->to_inevelocity(_dpdt);
			// 転がり速度，滑り速度の算出．（慣性座標系）
			//this->get_us_ur2(p, dpdt, us, ur);
			//this->get_us_ur3(p, er, us, ur);
			//Vector3d ey = Vector3d(0, 1, 0);	// 溝直交座標系のY方向単位ベクトル
			//Vector3d thXZ = this->RG->ine_to_XZcoord(this->BL->x);
			//Matrix3d R_th = this->RG->get_xyz2XYZ(thXZ[0]);
			//Matrix3d R_th_ = R_th.inverse();
			//Vector3d ey_ = R_th_ * ey;
			//Vector3d ey__ = this->RG->to_inevector(ey_);
			//double dot = ey_.dot(er);
			//Vector3d ex = Vector3d(1, 0, 0);
			//Vector3d e_th = this->RG->to_inevector(ex.cross(er));
			//Vector3d e_th_ = e_th / e_th.norm();
			//Vector3d ur_new = ey__ * (ey__.dot(ur));
			//double dot2 = ey__.dot(this->BL->x - this->RG->x);
			
			this->get_us_ur(p, us, ur);
			double ur_norm = ur.norm();

			// ボールにかかる滑り摩擦抗力の算出．
			double h = this->FT->calc(F_norm, this->E, Rx, Ry, this->lb_alpha, this->lb_eta, ur_norm, this->lm, b);
			h *= Tribology::ErtelGrubin(this->lb_eta, this->lb_beta, this->lb_k, ur_norm);
			lambda = h / this->sigma;
			fratio = Tribology::ForceRatio(lambda);
			this->calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);


			// ボールにかかる転がり摩擦抗力の算出．（慣性座標系）
			double Trf = this->RR->calc(Rx, Ry, this->BL->D, a, b, F_norm, this->E, ur_norm, fratio, this->lb_eta, this->lb_alpha, this->lb_beta, this->lb_k, this->lm);
			double Trh = this->HY->calc(this->fh, b, F_norm);
			Tbr = (Trf + Trh) * this->BL->calc_TorqueDirection(p, -ur);

			Fb = Fn + Fbs;
			Tb = Tbr + Tbs;
			Ti = this->RG->calc_Torque(p, -Fn) + Tis - Tbr;
		}
		else {
			// スライス片の結果をすべて0にする
			this->init_Sliceparam(i);
		}

		// 計算結果をメンバ変数に保存する．
		this->save(c, i, cos_alp, sin_alp, dx, a, b, p, Fn, Pmax, us, ur, Tbr, lambda, fratio, Fbs, Tbs, Fb, Tb, Ti);

		// 溝0と溝1の合計を出力する．
		Fbi += Fb;		// Fbi : ボール(b)がリング(i)から受ける力．
		Tbi += Tb;		// Tbi : ボール(b)がリング(i)から受けるモーメント．
		Fib += -Fb;		// Fib : リング(i)がボール(b)から受ける力．
		Tib += Ti;		// Tib : リング(i)がボール(b)から受けるモーメント．
	}
	return;
}

// 計算結果をメンバ変数から出力用変数に渡す
void B4P_BallRingPair::write(
	B4P_Out::BallRingPair&BRP	// out: (構造体): 玉-内外輪接触計算出力
) {

	this->thXZ = this->RG->ine_to_XZcoord(this->BL->x);
	BRP.th = this->thXZ[0];
	BRP.X = this->thXZ[1];
	BRP.Z = this->thXZ[2];
	Matrix3d xyz2XYZ = this->RG->get_xyz2XYZ(BRP.th);

	for (int ig = 0; ig < 2; ig++) {
		Vector3d Fr = this->GV[ig].Tr.norm() * this->BL->r_inv * -this->GV[ig].ur.normalized();
		for (int j = 0; j < 3; j++) {
			BRP.GV[ig].p[j] = this->GV[ig].p[j];
			BRP.GV[ig].Fn[j] = this->GV[ig].Fn[j];
			BRP.GV[ig].Fs[j] = this->GV[ig].Fs[j];
			BRP.GV[ig].Ts[j] = this->GV[ig].Ts[j];
			BRP.GV[ig].Fr[j] = Fr[j];
			BRP.GV[ig].Tr[j] = this->GV[ig].Tr[j];
			BRP.GV[ig].us[j] = this->GV[ig].us[j];
			BRP.GV[ig].ur[j] = this->GV[ig].ur[j];
		}
		Vector3d Fn_ = this->RG->ine_to_XZvector(this->GV[ig].Fn, xyz2XYZ);
		Vector3d Fs_ = this->RG->ine_to_XZvector(this->GV[ig].Fs, xyz2XYZ);
		Vector3d Ts_ = this->RG->ine_to_XZvector(this->GV[ig].Ts, xyz2XYZ);
		Vector3d Fr_ = this->RG->ine_to_XZvector(Fr, xyz2XYZ);
		Vector3d Tr_ = this->RG->ine_to_XZvector(this->GV[ig].Tr, xyz2XYZ);
		Vector3d us_ = this->RG->ine_to_XZvector(this->GV[ig].us, xyz2XYZ);
		Vector3d ur_ = this->RG->ine_to_XZvector(this->GV[ig].ur, xyz2XYZ);
		for (int j = 0; j < 3; j++) {
			BRP.GV[ig].Fn_[j] = Fn_[j];
			BRP.GV[ig].Fs_[j] = Fs_[j];
			BRP.GV[ig].Ts_[j] = Ts_[j];
			BRP.GV[ig].Fr_[j] = Fr_[j];
			BRP.GV[ig].Tr_[j] = Tr_[j];
			BRP.GV[ig].us_[j] = us_[j];
			BRP.GV[ig].ur_[j] = ur_[j];
		}
		Vector3d p_ = this->RG->ine_to_XZcoord(this->GV[ig].p);
		for (int j = 0; j < 2; j++)
			BRP.GV[ig].p_[j] = p_[j + 1];

		BRP.GV[ig].lambda = this->GV[ig].lambda;
		BRP.GV[ig].fratio = this->GV[ig].fratio;
		BRP.GV[ig].phi = this->ContactAngle(this->GV[ig].cos_alp, this->GV[ig].sin_alp);
		BRP.GV[ig].a = this->GV[ig].a;
		BRP.GV[ig].b = this->GV[ig].b;
		BRP.GV[ig].dx = this->GV[ig].dx;
		BRP.GV[ig].Pmax = this->GV[ig].Pmax;

		this->write_slice(ig, xyz2XYZ, BRP);
	}

	return;
}

// スライス片の結果を出力用構造体に格納
void B4P_BallRingPair::write_slice(
	int ig,						// in:	[-]: 溝番号
	Matrix3d xyz2XYZ,			// in:	[-]: 溝直交座標系変換行列
	B4P_Out::BallRingPair&BRP	// out:    : 玉-内外輪接触計算出力
) {
	for (int j = 0; j < this->GV[ig].msmax; j++) {
		BRP.GV[ig].SL[j].f_arr = this->GV[ig].SL[j].f_arr;
		BRP.GV[ig].SL[j].mu_cl = this->GV[ig].SL[j].mu_cl;
		BRP.GV[ig].SL[j].mu_tr = this->GV[ig].SL[j].mu_tr;
		Vector3d fs_ = this->RG->ine_to_XZvector(this->GV[ig].SL[j].fs, xyz2XYZ);
		Vector3d ts_ = this->RG->ine_to_XZvector(this->GV[ig].SL[j].ts, xyz2XYZ);
		Vector3d us_ = this->RG->ine_to_XZvector(this->GV[ig].SL[j].us, xyz2XYZ);
		Vector3d p_ = this->RG->ine_to_XZvector(this->GV[ig].SL[j].ps, xyz2XYZ);

		for (int k = 0; k < 3; k++) {
			BRP.GV[ig].SL[j].fs_[k] = fs_[k];
			BRP.GV[ig].SL[j].ts_[k] = ts_[k];
			BRP.GV[ig].SL[j].us_[k] = us_[k];
			BRP.GV[ig].SL[j].ps_[k] = p_[k];
		}
	}
	return;
}
double B4P_BallOuterRingPair::ContactAngle(double cos_alp, double sin_alp) {
	double alp = atan2(sin_alp, cos_alp);
	return alp;
}

double B4P_BallInnerRingPair::ContactAngle(double cos_alp, double sin_alp) {
	double alp = atan2(sin_alp, -cos_alp);
	return alp;
}

// 転がり速度を算出するメソッド（慣性座標系）
Vector3d B4P_BallRingPair::get_ur(
	const Vector3d&p		// 接触点位置．（慣性座標系）
) {
	Vector3d rwb = this->BL->w.cross(p - this->BL->x);
	Vector3d rwr = this->RG->w.cross(p - this->RG->x);
	Vector3d ur = 0.5 * (3 * this->RG->v - this->BL->v + rwb + rwr);
	Vector3d ur_ = ur;// this->BL->remove_Normal(p, ur);
	return ur_;
}

// 滑り速度を算出するメソッド（慣性座標系）
Vector3d B4P_BallRingPair::get_us(
	const Vector3d&p		// 接触点位置．（慣性座標系）
) {
	Vector3d ub = this->BL->surface_velocity(p);
	Vector3d ui = this->RG->surface_velocity(p);
	Vector3d us = ub - ui;
	Vector3d us_ = us;// this->BL->remove_Normal(p, us);
	return us_;
}

// 滑り速度と転がり速度を算出するメソッド（慣性座標系）．共通の変数が多いため高速化のために一体とした．
void B4P_BallRingPair::get_us_ur(
	const Vector3d&p,		// 接触点位置．（慣性座標系）
	Vector3d & us,			// out:	[-]:	滑り速度．（慣性座標系）
	Vector3d & ur			// out:	[-]:	転がり速度．（慣性座標系）
) {
	Vector3d dx = p - this->BL->x;
	Vector3d e = dx.normalized();
	Vector3d rwb = this->BL->w.cross(dx);
	Vector3d rwr = this->RG->w.cross(p - this->RG->x);
	Vector3d ur_ = 0.5 * (3 * this->RG->v - this->BL->v + rwb + rwr);
	ur = ur_;// -e//* e.dot(ur_);

	Vector3d ub = this->BL->v + rwb;
	Vector3d ui = this->RG->v + rwr;
	Vector3d us_ = ub - ui;
	us = us_;// -e; //* e.dot(us_);

	return;
}



// 転がり速度を算出するメソッド（慣性座標系）(第2案，接触点速度をより正確に計算しているが，演算量が増加)
void B4P_BallRingPair::get_us_ur2(
	const Vector3d&p,		// in:	[-]:	接触点位置．（慣性座標系）
	const Vector3d&er,		// in:	[-]:	リング中心から見たボール位相方向ベクトル．x成分=0．
	Vector3d & us,			// out:	[-]:	滑り速度．（慣性座標系）
	Vector3d & ur			// out:	[-]:	転がり速度．（慣性座標系）
) {
	Vector3d dx = p - this->BL->x;
	Vector3d e = dx.normalized();
	Vector3d rwb = this->BL->w.cross(dx);
	Vector3d rwr = this->RG->w.cross(p - this->RG->x);
	Vector3d bl_v = this->BL->v - this->RG->v;
	Vector3d e_a = this->RG->get_ax();
	Vector3d e_th = er.cross(e_a);				// 周方向ベクトル（径方向が12時の向きのとき3時の向き）
	Vector3d bl_vth = bl_v.dot(e_th)*e_th;

	Vector3d RB = this->BL->x - this->RG->x;
	// (軸方向速度にも処理をしなければならないが，まだ考慮できていない)
	Vector3d p_v = this->BL->v - this->RG->v + (bl_vth * dx.dot(er) / RB.dot(er));


	Vector3d ur_ = 0.5 * (this->BL->v + this->RG->v - 2 * p_v + rwb + rwr);
	ur = ur_;// -e * e.dot(ur_);

	Vector3d ub = this->BL->v + rwb;
	Vector3d ui = this->RG->v + rwr;
	Vector3d us_ = ub - ui;
	us = us_;// -e * e.dot(us_);

	return;
}

// B4P_In から読み込んで初期設定．
void B4P_BallOuterRingPair::init(const B4P_In&FI) {

	this->mu = FI.BOP.mu;		// 摩擦係数
	this->zeta = FI.BOP.dzeta;		// 減衰比

	this->init_Lubrication(FI.LB);
	this->init_Slice(FI.msmax);
	this->init_Tribology(FI.TB);

	return;
}

// B4P_In から読み込んで初期設定．
void B4P_BallInnerRingPair::init(const B4P_In&FI) {

	this->mu = FI.BIP.mu;		//摩擦係数
	this->zeta = FI.BIP.dzeta;		//減衰比

	this->init_Lubrication(FI.LB);
	this->init_Slice(FI.msmax);
	this->init_Tribology(FI.TB);

	return;
}

// 入力によって理論式をセットする．
void B4P_BallRingPair::init_Tribology(const B4P_In::Tribology&FI) {

	this->m = Tribology::ReducedMass(this->BL->m, this->RG->m);									// 換算質量
	this->E = Tribology::ReducedYoung(this->BL->E, this->BL->nu, this->RG->E, this->RG->nu);	// 等価ヤング率
	this->sigma = Tribology::CompositeRoughness(this->BL->sigmap, this->RG->sigmap);			// 合成粗さ

	switch (FI.rollingresistance) {
	case B4P_In::Tribology::RollingResistance::RollingResistanceNothing:
		this->RR = new Tribology::RollingResistanceNothing();
		break;
	case B4P_In::Tribology::RollingResistance::Aihara:
		this->RR = new Tribology::AiharaR();
		break;
	case B4P_In::Tribology::RollingResistance::Fujiwara:
		this->RR = new Tribology::Fujiwara();
		break;
	case B4P_In::Tribology::RollingResistance::Houpert:
		this->RR = new Tribology::Houpert();
		break;
	case B4P_In::Tribology::RollingResistance::GoksemAihara:
		this->RR = new Tribology::GoksemAihara();
		break;
	}

	this->HZ = new Tribology::BrewHamrock();

	switch (FI.coulomb) {
	case B4P_In::Tribology::CoulombNothing:
		this->CL = new Tribology::CoulombNothing();
		break;
	case B4P_In::Tribology::Tangent:
		this->CL = new Tribology::Tangent();
		this->clmb_s = FI.coulomb_slope;
		break;
	}

	switch (FI.filmThickness) {
	case B4P_In::Tribology::FilmThicknessNothing:
		this->FT = new Tribology::FilmThicknessNothing();
		break;
	case B4P_In::Tribology::HamrockDowsonHc:
		this->FT = new Tribology::HamrockDowsonHc();
		break;
	case B4P_In::Tribology::HamrockDowsonHmin:
		this->FT = new Tribology::HamrockDowsonHmin();
		break;
	}

	this->TR = new Tribology::AiharaT();

	this->DF = new Tribology::Tsuji();

	switch(FI.hysteresis) {
	case B4P_In::Tribology::Kakuta:
		this->HY = new Tribology::Kakuta();
		this->fh = FI.hysteresis_factor;
		break;
	case B4P_In::Tribology::HysteresisNothing:
		this->HY = new Tribology::HysteresisNothing();
		this->fh = 0;
		break;
	}
	
	return;
}

// 入力によって潤滑条件をセットする．
void B4P_BallRingPair::init_Lubrication(const B4P_In::Lubrication&FI) {

	this->lb_eta = FI.eta0;		// 粘度（Pa*s）
	this->lb_beta = FI.beta0;	// 温度粘度係数
	this->lb_k = FI.k0;			// 油熱伝導率
	this->lb_alpha = FI.alpha0;	// 圧力粘度係数
	this->lm = FI.lm0;			// メニスカス長さ

	return;
}

// 入力によって接触楕円分の変数を確保しておく．
void B4P_BallRingPair::init_Slice(int msmax) {

	for (int ig = 0; ig < 2; ig++) {
		int n = msmax;
		this->GV[ig].msmax = n;
		this->GV[ig].ps = new Vector3d[n];
		this->GV[ig].ratio_slice = new double[n];
		this->GV[ig].SL = new Slice[n];
		// （必要ないと思うが，）スライス片の数値を初期化
		this->init_Sliceparam(ig);
		// Hertz接触理論からスライス荷重を求める．
		Tribology::SliceForceRatio(n, this->GV[ig].ratio_slice);
	}
	return;
}

B4P_BallRingPair::B4P_BallRingPair() {
	for (int i = 0; i < 2; i++) {
		this->GV[i].ps = NULL;
		this->GV[i].ratio_slice = NULL;
	}
	return;
}

B4P_BallRingPair::~B4P_BallRingPair() {
	for (int i = 0; i < 2; i++) {
		if (this->GV[i].ps != NULL)
			delete[] this->GV[i].ps;
		if (this->GV[i].ratio_slice != NULL)
			delete[] this->GV[i].ratio_slice;
		if (this->GV[i].SL != NULL)
			delete[] this->GV[i].SL;
	}
	return;
}

// STEP1 静解析　剛性計算：玉座標の代入
void B4P_BallOuterRingPair::set_eta_stf(const Vector3d & eta) {
	Vector3d x = this->RG->XZ_to_inecoord(eta);
	this->BL->x = x;
	return;
}


// STEP1 静解析　剛性計算：玉ー内外輪間の接触荷重の計算
void B4P_BallRingPair::get_Fstf(
	Vector3d&Fbi,	// out: [N]:	転動体が内輪（or外輪）から受ける力 
	Vector3d&Fib,	// out: [N]:	内輪（or外輪）が転動体から受ける力
	Vector3d&Tib	// out:	[Nm]:	内輪（or外輪）が転動体から受けるトルク
) {
	Fbi = Vector3d::Zero();
	Fib = Tib = Vector3d::Zero();

	Vector3d bl_x = this->RG->to_mycoord(this->BL->x);
	Vector3d bl_v = this->RG->to_myvelocity(this->BL->v);
	Vector3d Tbi;

	// groove と ball の接触計算．
	for (int i = 0; i < 2; i++) {

		// 必要な戻り値の確保．
		double dx, cos_alp, sin_alp, a, b, Pmax, lambda, fratio;
		dx = a = b = Pmax = lambda = fratio = 0;
		Vector3d er, eg, p, Fn, ur, us, Tr, Fs, Ts, Fi, Fb, Tb, Ti;
		er = eg = p = Fn = ur = us = Tr = Fs = Ts = Fi = Fb = Tb = Ti = Vector3d::Zero();

		// 接していているか判定．
		bool c = this->how_Contact(i, bl_x, er, eg, dx);
		if (c) {
			double Rx, Ry, k;
			this->calc_Hertz(i, bl_x, er, eg, dx, Rx, Ry, p, cos_alp, sin_alp, a, b, k);
			double F_norm = k * pow(dx, 1.5);
			Fb = F_norm * this->RG->to_inevector(-eg);
			Fi = F_norm * this->RG->to_inevector(eg);
			Pmax = 1.5 * F_norm / (Numeric::pi * a * b);
			Ti = this->RG->calc_Torque(p, Fi);	// 玉にかかる力の反作用から溝にかかるトルクを計算する．
		}
		// 接触していない場合でも弱いばねで接続
		else {
			double k = 1e-10;
			double F_norm = k * dx;
			Fb = F_norm * this->RG->to_inevector(-eg);
			Fi = F_norm * this->RG->to_inevector(eg);
			Ti = this->RG->calc_Torque(p, Fi);	// 玉にかかる力の反作用から溝にかかるトルクを計算する．
		}
		// 計算結果をメンバ変数に保存する．
		//this->save(c, i, cos_alp, sin_alp, a, b, p, Fn, Pmax, us, ur, Tr, lambda, fratio, Fs, Ts, Fb, Tb, Ti);

		// 溝0と溝1の合計を出力する．
		Fbi += Fb;		// Fbi : ボール(b)がリング(i)から受ける力．
		Fib += Fi;		// Fib : リング(i)がボール(b)から受ける力．
		Tib += Ti;
	}

	return;
}


// 各スライスの荷重を出力
//VectorXd B4P_BallRingPair::make_SliceParam(void) {
//	int NParam = 15;
//	int Nslice = this->GV[0].msmax + this->GV[1].msmax;
//	VectorXd Param(Nslice * NParam);
//	for (int i = 0; i < 2; i++) {
//		int Ni = i * NParam * this->GV[0].msmax;
//		for (int k = 0; k < this->GV[i].msmax; k++) {
//			int Nik = Ni + NParam * k;
//			for (int j = 0; j < 3; j++) {
//				int Nikj = Nik + j;
//				Param(Nikj + 0) = this->GV[i].ps[k][j];
//				Param(Nikj + 3) = this->GV[i].Fs_slice[k][j];
//				Param(Nikj + 6) = this->GV[i].Ts_slice[k][j];
//				Param(Nikj + 9) = this->GV[i].us_slice[k][j];
//			}
//			Param(Nik + 12) = this->GV[i].fn_slice[k];
//			Param(Nik + 13) = this->GV[i].mutr_slice[k];
//			Param(Nik + 14) = this->GV[i].mucl_slice[k];
//		}
//	}
//	return Param;
//}



//Vector3d ub = this->BL->surface_velocity(p);
//Vector3d ui = this->RG->surface_velocity(p);
//Vector3d ur = 0.5 * (ub + ui);


// 部材の材料と接触部剛性から，ダンピング係数を求める．"15_非線形バネと衝突のモデル化"参考．
//// 計算式は"日本惑星科学会誌, 2004, 13.4: 233-240."の非線形ばねへの手計算による適応． 
//double B4P_BallRingPair::get_c
//(				// out:	[Ns/m]:	ダンピング係数．
//	double k	// in:	[N/m]:	接触部剛性．
//) {
//	double c = 2 * this->zeta * sqrt(1.5 * this->m * k);
//	return c;
//}




//// 遠心力を算出するメソッド（静解析用）戻り値：慣性座標系[N]
//Vector3d B4P_BallRingPair::get_centf() {
//
//	// 内外輪座標系(＝内外輪中心基準)の玉の位置ベクトル
//	Vector3d r_rg = this->RG->to_mycoord(this->BL->x);
//
//	// x成分は関係ないので0に設定
//	r_rg.x() = 0;
//
//	// r_のノルムと単位ベクトルを定義
//	double rrg_nor = r_rg.norm();
//	Vector3d e_rg = r_rg / rrg_nor;
//
//	// 内外輪座標系の玉速度の接線方向成分を算出
//	Vector3d v_rg = this->RG->to_mycoord(this->BL->v);
//	v_rg.x() = 0;
//	double vrg_th = (v_rg.cross(e_rg)).norm();
//
//	// 遠心力の計算式(F=m*v^2/r)を適用し，内外輪座標系上の遠心力を出力
//	Vector3d F_rg = e_rg * this->BL->m * vrg_th * vrg_th / rrg_nor;
//
//	// 内外輪座標系の遠心力を慣性座標系に変換
//	Vector3d F_inr = this->RG->to_inevector(F_rg);
//
//	// 慣性座標系における遠心力を出力
//	return F_inr;
//}


//// 転がり速度を算出するメソッド（慣性座標系）(第xx案，多分間違っている)
//void B4P_BallRingPair::get_us_ur2(
//	const Vector3d&p,		// in:	[-]:	接触点位置．（慣性座標系）
//	const Vector3d&dpdt,		// in:	[-]:	接触点位置．（慣性座標系）
//	Vector3d & us,			// out:	[-]:	滑り速度．（慣性座標系）
//	Vector3d & ur			// out:	[-]:	転がり速度．（慣性座標系）
//) {
//	Vector3d dx = p - this->BL->x;
//	Vector3d e = dx.normalized();
//	Vector3d rwb = this->BL->w.cross(dx);
//	Vector3d rwr = this->RG->w.cross(p - this->RG->x);
//	Vector3d bl_v = this->BL->v - this->RG->v;
//	//Vector3d e_a = this->RG->get_ax();
//
//	Vector3d RB = this->BL->x - this->RG->x;
//	
//
//
//	Vector3d ur_ = 0.5 * (this->BL->v + this->RG->v - 2 * dpdt + rwb + rwr);
//	ur = ur_-e * e.dot(ur_);
//
//	Vector3d ub = this->BL->v + rwb;
//	Vector3d ui = this->RG->v + rwr;
//	Vector3d us_ = ub - ui;
//	us = us_;// -e //* e.dot(us_);
//
//	return;
//}