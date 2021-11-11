#include "B4P_BallCagePair.h"

// 玉・保持器オブジェクトをメンバ変数に代入
void B4P_BallCagePair::link(Ball*BL, bal_Cage*CG, int np) {

	this->BL = BL;
	this->CG = CG;
	this->np = np;

	return;
}

// パラメータをメンバ変数に代入
void B4P_BallCagePair::init(const B4P_In&FI) {

	this->mu = FI.BCP.mu;
	this->zeta = FI.BCP.dzeta;
	if (false) {
		this->DF = new Tribology::Tsuji();
		Tribology::BrewHamrock BH;
		double Rx =  1. / (1. / FI.Snap.R + 1. / this->BL->r);
		double E = 1. / (1. / FI.Cage.E + 1. / FI.BL[0].E);
		double k, a, b;
		BH.calc(Rx, Rx, 1, E, k, a, b);
		printf("%f\n", k);
	}
	else {
		this->DF = new Tribology::KelvinVoigt();
	}
	this->m = Tribology::ReducedMass(this->BL->m, this->CG->m);

	return;
}

// 玉-保持器間に働く荷重・トルクを計算
void B4P_BallCagePair::calc_force(Vector3d&Fbc, Vector3d&Tbc, Vector3d&Fcb, Vector3d&Tcb) {

	// 戻り値の宣言および初期化．
	double k[_MAX_CONTACT_], dx[_MAX_CONTACT_]; 
	int ptt[_MAX_CONTACT_];	// 接触パターン
	Vector3d p[_MAX_CONTACT_], Fn[_MAX_CONTACT_], Fs[_MAX_CONTACT_];	// 接触点位置における剛性．
	for (int i = 0; i < _MAX_CONTACT_; i++) {
		p[i] = Fn[i] = Fs[i] = Vector3d::Zero();
		ptt[i] = 0;
	}
	// 接触点位置を保持器に判別してもらう．
	int num_cont = this->CG->get_ContactPoint(this->BL->x, this->BL->r, this->np, p, k, ptt);

	// ループで積算を求めていくが，その前に初期化を行う．
	Fbc = Tbc = Fcb = Tcb = Vector3d::Zero();

	// 接触点の数だけ荷重・摩擦計算をし，加算していく．
	for (int i = 0; i < num_cont; i++) {

		// 接触荷重の計算
		double _dx = 0, dv = 0;
		Vector3d edir = Vector3d::Zero();
		Vector3d v = this->CG->surface_velocity(p[i]);
		bool c = this->BL->calc_Contact(p[i], v, _dx ,dv, edir);

		// 接触しない場合は後の処理を省略．
		if (!c) 
			continue;

		// 減衰も考慮した垂直荷重の計算
		double fn = this->DF->calc(k[i], this->zeta, this->m, dv, _dx);
		Vector3d _Fn = fn * -edir;

		// 滑り摩擦の計算
		Vector3d us = this->get_us(p[i]);
		double us_norm = us.norm();
		Vector3d _Fs = Vector3d::Zero();
		if (us_norm != 0.0) 
			_Fs = this->mu * _Fn.norm() * -us / us_norm;
		
		Fn[i] = _Fn;
		Fs[i] = _Fs;
		dx[i] = _dx;
		Fbc += _Fn + _Fs;
		Tbc += this->BL->calc_Torque(p[i], Fbc);
		Tcb += this->CG->calc_Torque(p[i], -Fbc);
	}
	// 計算結果をメンバ変数に格納
	for (int i = 0; i < _MAX_CONTACT_; i++) {
		this->p[i] = p[i];
		this->Fn[i] = Fn[i];
		this->Fs[i] = Fs[i];
		this->dx[i] = dx[i];
		this->ptt[i] = ptt[i];
	}
	Fcb = -Fbc;

	return;
}


// 滑り速度を算出するメソッド（慣性座標系）
Vector3d B4P_BallCagePair::get_us(
	const Vector3d&p		// 接触点位置．（慣性座標系）
) {
	Vector3d ub = this->BL->surface_velocity(p);
	Vector3d uc = this->CG->surface_velocity(p);
	Vector3d us = ub - uc;
	Vector3d us_ = this->BL->remove_Normal(p, us);
	return us_;
}

// 最新の状態の，玉における，接触点位置，接触荷重ベクトル，滑り摩擦力ベクトル，転がり摩擦力ベクトル，滑り速度ベクトル，転がり速度ベクトル，ラムダ値，荷重支持割合，接触角の余弦・正弦，接触楕円長半径・短半径
void B4P_BallCagePair::save(B4P_Out::BallCagePair & BCP) {

	for (int i = 0; i < _MAX_CONTACT_; i++) {
		Vector3d p_ = this->CG->to_mycoord(this->p[i]);
		Vector3d Fn_ = this->CG->to_myvector(this->Fn[i]);
		Vector3d Fs_ = this->CG->to_myvector(this->Fs[i]);

		for (int j = 0; j < 3; j++) {
			BCP.CP[i].p[j] = this->p[i][j];
			BCP.CP[i].Fn[j] = this->Fn[i][j];
			BCP.CP[i].Fs[j] = this->Fs[i][j];
			BCP.CP[i].p_[j] = p_[j];
			BCP.CP[i].Fn_[j] = Fn_[j];
			BCP.CP[i].Fs_[j] = Fs_[j];
		}
		BCP.CP[i].dx = this->dx[i];
		BCP.CP[i].ptt = this->ptt[i];
	}
	return;
}

