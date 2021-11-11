#include "bal_Cage.h"

// 抽象保持器クラス．
void bal_Cage::Nothing(void) {
}

// 冠型保持器クラス．
void bal_SnapCage::init(const B4P_In&FI) {

	this->ro = FI.Snap.dout * 0.5;
	this->ri = FI.Snap.din * 0.5;
	this->m = FI.Cage.m;
	this->pcd = this->ro + this->ri + FI.Snap.jc;
	this->m_inv = 1.0 / this->m;
	for (int i = 0; i < 3; i++)
		this->xc[i] = FI.Cage.rmg0[i];
	this->Z = FI.ballnum;
	this->I = Vector3d(FI.Cage.Ix, FI.Cage.Iyz, FI.Cage.Iyz);
	this->I_inv = Vector3d(1 / FI.Cage.Ix, 1 / FI.Cage.Iyz, 1 / FI.Cage.Iyz);

	double x; ArrayXd th(this->Z + 1), y(this->Z + 1), z(this->Z + 1);
	th = ArrayXd::LinSpaced(this->Z + 1, 0, 2 * Numeric::pi);
	x = 0;
	y = 0.5 * this->pcd *-th.sin();
	z = 0.5 * this->pcd * th.cos();

	// ポケット配列を動的確保．
	this->PK = new Pocket[FI.ballnum];

	// 全ポケットを初期化．
	for (int i = 0; i < Z; i++) {
		Vector3d ax = this->ax0;		// 軸方向
		this->PK[i].R = FI.Snap.R;
		this->PK[i].x = Vector3d(x, y[i], z[i]);
		this->PK[i].er = this->PK[i].x.normalized();
		this->PK[i].eth = this->PK[i].er.cross(ax);
		this->PK[i].x_norm = this->PK[i].x.norm();

		this->PK[i].ropen = FI.Snap.ropen;
		this->PK[i].kface = FI.Snap.Kface;			// 面接触剛性[N/m]
		this->PK[i].kci = FI.Snap.Kcornerin;		// 角接触剛性[N/m]
		this->PK[i].kco = FI.Snap.Kcornerout;		// 角接触剛性[N/m]
		this->PK[i].kopen = FI.Snap.Kopen;			// 開口部剛性[N/m]
		this->PK[i].kedgei = FI.Snap.Kedgein;			// エッジ剛性[N/m]	
		this->PK[i].kedgeo = FI.Snap.Kedgeout;			// エッジ剛性[N/m]

		// 開口部開き角α
		this->PK[i].cos_alp = sqrt(1.0 - FI.Snap.ropen*FI.Snap.ropen / this->PK[i].R / this->PK[i].R);
		this->PK[i].h0 = this->PK[i].R * this->PK[i].cos_alp;

		// 開口部開き角γ
		double cosi = Numeric::Law_of_cosines(this->PK[i].ropen, ri, this->PK[i].x_norm);
		double coso = Numeric::Law_of_cosines(this->PK[i].ropen, ro, this->PK[i].x_norm);
		double sini = sqrt(1.0 - cosi * cosi);
		double sino = sqrt(1.0 - coso * coso);
		this->PK[i].cos_gammai = (ri * cosi - this->PK[i].x_norm) / sqrt(pow(ri * cosi - this->PK[i].x_norm, 2) + pow(ri * sini, 2));
		this->PK[i].cos_gammao = (ro * coso - this->PK[i].x_norm) / sqrt(pow(ro * coso - this->PK[i].x_norm, 2) + pow(ro * sino, 2));

		// 垂直高さの計算
		double cosim = Numeric::Law_of_cosines(this->PK[i].R, ri, this->PK[i].x_norm);
		this->PK[i].zimax = this->ri;
		this->PK[i].zimin = this->ri * cosim;

		double cosom = Numeric::Law_of_cosines(this->PK[i].R, ro, this->PK[i].x_norm);
		this->PK[i].zomax = this->ro;
		this->PK[i].zomin = this->ro * cosom;

		this->PK[i].ziave = (this->PK[i].zimax + this->PK[i].zimin) * 0.5;
		this->PK[i].ziamp = (this->PK[i].zimax - this->PK[i].zimin) * 0.5;

		this->PK[i].zoave = (this->PK[i].zomax + this->PK[i].zomin) * 0.5;
		this->PK[i].zoamp = (this->PK[i].zomax - this->PK[i].zomin) * 0.5;

		// コーナーの座標の設定．
		this->PK[i].corni0 = this->PK[i].h0 * ax - ri * sini * this->PK[i].eth + ri * cosi * this->PK[i].er;
		this->PK[i].corni1 = this->PK[i].h0 * ax + ri * sini * this->PK[i].eth + ri * cosi * this->PK[i].er;
		this->PK[i].corno0 = this->PK[i].h0 * ax - ro * sino * this->PK[i].eth + ro * coso * this->PK[i].er;
		this->PK[i].corno1 = this->PK[i].h0 * ax + ro * sino * this->PK[i].eth + ro * coso * this->PK[i].er;


	}
	return;
}

// ポケット中心を基準とした玉方向ベクトルから接触パターンを判定（接触の有無は判定しない）
bal_SnapCage::ContactPattern bal_SnapCage::how_Contact
(							// out:	[-]:	接触パターン．
	const Vector3d&BL_x_,	// in:	[m]:	ボールの位置（保持器幾何中心座標系）．
	int np,					// in:	[-]:	そのボールの対応するポケット番号．
	Vector3d&edir,			// out: [-]:	ポケット中心から玉中心へ向かう単位ベクトル．
	double&zo_th,			// out: [m]:	保持器中心からポケット中心方向で見た，現在の方位角における保持器外エッジの高さ．
	double&zi_th, 			// out: [m]:	保持器中心からポケット中心方向で見た，現在の方位角における保持器内エッジの高さ．
	double&cos_th, 			// out: [-]:	軸方向を12時に見立てた際の方位角のcosine．
	double&sin_th, 			// out: [-]:	軸方向を12時に見立てた際の方位角のsine．
	double&cos_gamma, 		// out: [-]:	軸方向からみたポケットを中心とした玉の位相角のcosine．
	double&sin_gamma		// out: [-]:	軸方向からみたポケットを中心とした玉の位相角のsine．
) {
	// 前準備としてポケット中心を基準とした玉方向ベクトル edir を計算．（保持器幾何中心座標系）
	Vector3d pb = BL_x_ - this->PK[np].x;
	double pb_norm = pb.norm();

	// 0除算の場合は例外処理（ポケット中心に玉がある場合）
	if (pb_norm == 0.0)
		return exception;
	edir = pb / pb_norm;

	// (1) 玉方位角αを導出（軸方向を基準としたときの方位角，開口部当たり判定に使用）
	// 方位角αの余弦を計算．
	double   cos_alp = this->ax0.dot(edir);	// 注：α=[0,180]なのでcosで一意に決まる．

	// (2) 玉方位角γを導出（"軸方向からみた"径方向を基準としたときの方位角，角当たり判定に使用）
	// 外積を用いて径方向と γ+90° をなすベクトルを計算
	Vector3d edirax = edir.cross(this->ax0);
	double edirax_norm = edirax.norm();

	// 0除算の場合は例外処理（軸とポケット中心-ボール中心が並行）
	if (edirax_norm == 0.0)
		return ContactPattern::exception;

	Vector3d azi = edirax / edirax_norm;
	cos_gamma = this->PK[np].eth.dot(azi);
	sin_gamma = -this->PK[np].er.dot(azi);

	// (3) 垂直高さの計算（エッジ当たりに使用）
	// 保持器幾何中心・ポケット中心・ボール中心の張る三角形の法線ベクトル．（保持器幾何中心座標系）
	Vector3d eredir = this->PK[np].er.cross(edir);
	double   eredir_norm = eredir.norm();

	// 0除算の場合は例外処理（保持器幾何中心・ポケット中心・ボール中心が一列）
	if (eredir_norm == 0.0)
		return exception;

	Vector3d en = eredir / eredir_norm;

	// 方位角 θ の sine.（12時が保持器幾何中心座標系の[1,0,0]方向，3時が周方向．）
	sin_th = -this->ax0.dot(en);
	double cos_2th = 1.0 - 2 * sin_th * sin_th;
	cos_th = this->PK[np].eth.dot(en);

	// 近似式によって方位による垂直高さを仮定．（玉径に対して誤差 0.01% 程度を確認済み．）
	zi_th = this->PK[np].ziave + this->PK[np].ziamp * cos_2th;
	zo_th = this->PK[np].zoave + this->PK[np].zoamp * cos_2th;

	// 接触判定候補の垂直高さ L の算出．
	double cos_beta = this->PK[np].er.dot(edir);
	double L = this->PK[np].x_norm + this->PK[np].R * cos_beta;

	// 最後に場合分けを実施．
	if (cos_alp < this->PK[np].cos_alp) {
		if (L < zi_th)
			return edgei;
		else if (L > zo_th)
			return edgeo;
		else
			return face;
	}
	if (cos_gamma < this->PK[np].cos_gammai)
		return corneri;
	if (cos_gamma > this->PK[np].cos_gammao)
		return cornero;

	return aperture;
}

// 接触パターンに応じて接触点候補を計算するメソッド．保持器は剛体と仮定．（保持器幾何中心座標系）
int bal_SnapCage::where_Contact
(							// out:	[-]:	接触点数．
	ContactPattern c,		// in:	[-]:	接触パターン．
	int np,					// in:	[-]:	そのボールの対応するポケット番号．
	const Vector3d&dir, 	// in:	[-]:	ポケット中心から玉中心へ向かう単位ベクトル．
	double zo_th, 			// in:	[m]:	保持器中心からポケット中心方向で見た，現在の方位角における保持器外エッジの高さ．
	double zi_th, 			// in:	[m]:	保持器中心からポケット中心方向で見た，現在の方位角における保持器内エッジの高さ．
	double cos_th, 			// in:	[-]:	軸方向を12時に見立てた際の方位角のcosine．
	double sin_th, 			// in:	[-]:	軸方向を12時に見立てた際の方位角のsine．
	double cos_gamma, 		// in:	[-]:	軸方向からみたポケットを中心とした玉の位相角のcosine．
	double sin_gamma, 		// in:	[-]:	軸方向からみたポケットを中心とした玉の位相角のsine．
	Vector3d * x_,			// out:	[m]:	接触点位置（保持器幾何中心座標系）．接触点数分だけポインタで返
	double * k	   			// out:	[N/m]:	各接触点における剛性．接触点数分だけポインタで返す．
) {
	switch (c) {
	case face: {
		x_[0] = this->PK[np].x + dir * this->PK[np].R;
		k[0] = this->PK[np].kface;
		return 1;
	}
	case edgeo: {
		double R_cs = sqrt(Numeric::Square(this->PK[np].R) - Numeric::Square(zo_th - this->PK[np].x_norm));
		x_[0] = this->ax0    * cos_th * R_cs
			+ this->PK[np].eth * sin_th * R_cs
			+ this->PK[np].er * zo_th;
		k[0] = this->PK[np].kedgeo;
		x_[1] = - this->ax0    * cos_th * R_cs
			- this->PK[np].eth * sin_th * R_cs
			+ this->PK[np].er * zo_th;
		k[1] = this->PK[np].kedgeo;
		return 2;
	}
	case edgei: {
		double R_cs = sqrt(Numeric::Square(this->PK[np].R) - Numeric::Square(zi_th - this->PK[np].x_norm));
		x_[0] = this->ax0    * cos_th * R_cs
			+ this->PK[np].eth * sin_th * R_cs
			+ this->PK[np].er * zi_th;
		k[0] = this->PK[np].kedgei;
		x_[1] = - this->ax0    * cos_th * R_cs
			- this->PK[np].eth * sin_th * R_cs
			+ this->PK[np].er * zi_th;
		k[1] = this->PK[np].kedgei;
		return 2;
	}
	case aperture: {
		x_[0] = this->PK[np].x
			+ this->ax0      * this->PK[np].h0
			+ this->PK[np].eth * sin_gamma * this->PK[np].ropen
			+ this->PK[np].er * cos_gamma * this->PK[np].ropen;
		k[0] = this->PK[np].kopen;
		x_[1] = this->PK[np].x
			+ this->ax0      * this->PK[np].h0
			- this->PK[np].eth * sin_gamma * this->PK[np].ropen
			- this->PK[np].er * cos_gamma * this->PK[np].ropen;
		k[1] = this->PK[np].kopen;
		return 2;
	}
	case cornero: {
		x_[0] = this->PK[np].corno0;
		x_[1] = this->PK[np].corno1;
		k[0] = this->PK[np].kco;
		k[1] = this->PK[np].kco;
		return 2;
	}
	case corneri: {
		x_[0] = this->PK[np].corni0;
		x_[1] = this->PK[np].corni1;
		k[0] = this->PK[np].kci;
		k[1] = this->PK[np].kci;
		return 2;
	}
	case exception:
		return 0;

	default:
		return 0;
	}
}


// ボールの位置と半径によって接触点位置と剛性を返すメソッド．この情報をもとに外部で荷重摩擦計算を行う．
int bal_SnapCage::get_ContactPoint(
	const Vector3d&BL_x,	// in:	[m]:	ボールの位置（慣性座標系）．
	double         BL_r,	// in:	[m]:	ボールの半径． 
	int            np, 		// in:	[-]:	そのボールの対応するポケット番号．
	Vector3d      *x, 		// out:	[m]:	接触点位置（慣性座標系）．接触点数分だけポインタで返す．
	double        *k, 		// out:	[N/m]:	各接触点における剛性．接触点数分だけポインタで返す．
	int			  *ptt		// out: [-]:	接触パターン(int型に変換して出力)
) {
	// ボール位置を求める．（保持器幾何中心座標系）
	Vector3d BL_x_ = this->to_mycoord(BL_x) - this->xc;

	Vector3d dir = Vector3d::Zero();
	double zo_th, zi_th, cos_th, sin_th, cos_gamma, sin_gamma;
	zo_th = zi_th = cos_th = sin_th = cos_gamma = sin_gamma = 0.0;

	ContactPattern c = this->how_Contact(BL_x_, np, dir, zo_th, zi_th, cos_th, sin_th, cos_gamma, sin_gamma);

	Vector3d x_[_MAX_CONTACT_];	// 接触点位置の配列．（保持器幾何中心座標系）
	int nc = this->where_Contact(c, np, dir, zo_th, zi_th, cos_th, sin_th, cos_gamma, sin_gamma, x_, k);

	// 保持器幾何中心座標系 → 保持器座標系 → 慣性座標系に変換．
	for (int i = 0; i < nc; i++) {
		x[i] = this->to_inecoord(x_[i] + this->xc);
		ptt[i] = static_cast<int>(c);
	}
	return nc;
}

bal_SnapCage::bal_SnapCage() {
	this->PK = NULL;
	return;
}

bal_SnapCage::~bal_SnapCage() {
	if (this->PK != NULL)
		delete[] this->PK;
	return;
}



//冠型保持器重量自動計算
//double bal_SnapCage::calc_snapmass(double Snap_den_, double Snap_rout_, double Snap_rin_, double Snap_h_, double ballpcd_, double balldia_, double clr_, double ballnum_){
//	double rp3 = 0.5*ballpcd_+0.5*balldia_+0.25*clr_;//外輪半径
//	double mc;
//	/* snap cage 1.4465 is modification factor */
//	mc=1.4465*Snap_den_*Numeric::pi*((Snap_rout_*Snap_rout_-Snap_rin_*Snap_rin_)*Snap_h_
//		-(double)ballnum_*(Snap_rout_-Snap_rin_)*rp3*rp3
//		*(0.5+0.5*sin(2.0*asin((Snap_h_)/rp3))+asin((Snap_h_)/rp3)));
//	return mc;
//}

