/*******************************************************************************
!								"Tribology.cpp"
!													2018/11/27	[Core-T]	楠崎
!
!
!
!
!
!*******************************************************************************/

#include "Tribology.h"

// 接触荷重のスライス分割を積分により計算．
void Tribology::SliceForceRatio
(
	int msmax,			// in:	[-]:	接触楕円長径方向にスライスする数．
	double*ForceRatio	// out:	[N]:	接触荷重割合分布．配列サイズ＝msmax．
	) {
	for (int i = 0; i < msmax; i++) {
		double x1 = -1 + 2 * i / (double)msmax;					// i番目の要素のx座標を短冊の幅で割った数値
		double x2 = -1 + 2 * (i + 1) / (double)msmax;
		ForceRatio[i] = 0.25 * (3 * (x2 - x1) - (Numeric::Cube(x2) - Numeric::Cube(x1)));
	}
	return;
}

// 接触荷重のスライス分割を積分により計算．
void Tribology::SliceForceRatio2d
(
	int nr,		// in:	[-]:	接触楕円r方向にスライスする数．
	int nt,		// in:	[-]:	接触楕円θ方向にスライスする数．
	double*r	// out:	[N]:	接触荷重割合分布．配列サイズ:nr．
) {
	double dt = 2 * Numeric::pi / nt;

	for (int i = 0; i < nr; i++) {
		double x0 = double(i)     / nr;
		double r0 = -pow((1 - x0 * x0), 1.5) / 2 / Numeric::pi;

		double x1 = double(i + 1) / nr;
		double r1 = -pow((1 - x1 * x1), 1.5) / 2 / Numeric::pi;
		
		r[i] = (r1 - r0) * dt;
	}
	return;
}

// 等価ヤング率を算出する．
double Tribology::ReducedYoung
(						// out:	[Pa]:	等価ヤング率．
	double E0,			// in:	[Pa]:	物体0のヤング率．
	double m0,			// in:	[-]:	物体0のポアソン比．
	double E1,			// in:	[Pa]:	物体1のヤング率． 
	double m1			// in:	[-]:	物体0のポアソン比．
) {
	double E = 2.0 / ((1.0 - m0*m0)/E0 + (1.0 - m1*m1)/E1);
	return E;
}

// 換算質量を算出する．
double Tribology::ReducedMass
(						// out:	[kg]:	換算質量．
	double m0,			// in:	[kg]:	物体0の質量．
	double m1			// in:	[kg]:	物体1の質量．
) {
	double M = 1.0 / (1.0 / m0 + 1.0 / m1);
	return M;
}

// 合成粗さを算出する．
double Tribology::CompositeRoughness
(						// out:	[m]:	合成粗さ．
	double s0,			// in:	[m]:	物体0の粗さ．
	double s1			// in:	[m]:	物体1の粗さ．
) {
	double s = sqrt(s0 * s0 + s1 * s1);
	return s;
}

// ラムダ値によって油膜荷重割合を算出する．
double Tribology::ForceRatio
(						// out:	[-]:	油膜荷重支持割合．
	double lambda		// in:	[-]:	ラムダ値．
) {
	if (lambda < 0)
		return 0.0;

	return 1.0 - exp(-1.8 * pow(lambda, 1.2));
}
// 玉-溝接触時の接触楕円半径からころ有効長さ（の半分）を計算（研内第457号）
double Tribology::RollerLength(
	double a,	// 接触部長半径 [m]
	double Rx,	// 転がり方向等価半径 [m]
	double Ry,	// 転がり方向と直交する向きの等価半径 [m]
	double xio	// メニスカス長さ [m]
) {

	// メニスカス長さが未定義のときはaの1.2倍をころ長さとする．
	if (xio < 0) {
		double wi = 1.2 * a;
		return wi;
	}
	double Xio = xio / Rx;
	double A = a / Rx;
	double Hio = 0.5 * Xio * Xio;
	double ayw = Ry / Rx;

	double wi = Ry * sqrt(2 * ayw * Hio + A * A);
	// ころ有効長さがaの1.2倍を超えないようにする（理由は研内参照）
	if (wi > 1.2 * a) {
		wi = 1.2 * a;
	}

	return wi;
}

// タンジェントカーブの摩擦係数．
double Tribology::Tangent::calc(double mu0, double us, double ur, double s) {

	return mu0 * 2.0 / Numeric::pi * atan(s * us / ur);
}

// クーロン摩擦がない状態を定義．
double Tribology::CoulombNothing::calc(double mu0, double us, double ur, double s) {
	return 0.0;
}


double Tribology::SimpleCoulomb::calc(double mu0, double us, double ur, double s) {
	return mu0;
}

// Brew Hamrock の近似計算．教科書の値との一致を確認．
void Tribology::BrewHamrock::calc
(
	double Rx,			// in:	[m]:	周方向等価曲率．
	double Ry,			// in:	[m]:	転がり方向等価曲率．
	double dx,			// in:	[m]:	ボール食い込み量．
	double E,			// in:	[Pa]:	等価ヤング率．
	double&k,			// out:	[N/m^1.5]:	非線形剛性．
	double&a,			// out:	[m]:	接触楕円長径．
	double&b			// out:	[m]:	接触楕円短径．
) {
	double sum_rho = 1. / Rx  + 1. / Ry;
	double Kk = 1.5277 + 0.6023 * log(Ry/Rx);										// 第一種完全積分
	double Ek = 1.0003 + 0.5968 * Rx / Ry;											// 第二種完全積分
	double mu = 0.87959 * pow(Ek, Numeric::ot) * pow((Ry/Rx), 0.4240);				// 係数 μ

	double k_el = 1 / (1.0339 * pow(Ry/ Rx, 0.636));								// 楕円率 k_el 
	double nu = pow(2 * Ek * k_el / Numeric::pi, Numeric::ot);						// 係数 ν

	double K2pimu = 2 * Kk / Numeric::pi / mu;										// 係数 K/2πμ

	k = pow(pow(1.125 / E / E * sum_rho, Numeric::ot) * K2pimu, -1.5);				// 弾性係数の算出
	a = pow(6 * Ek * k / E / sum_rho / Numeric::pi / k_el / k_el , Numeric::ot) * sqrt(dx);	// 長半径 a の算出
	b = k_el * a;																	// 短半径 b の算出
	return;
}

// Hydrodynamic force and moment in pure rolling lubricated contacts. Part 2: Point contacts に基づき，トルク計算を行う．
double Tribology::Houpert::calc
(						// out:	[Nm]:	トルク．
	double Rx,			// in:	[m]:	転がり方向等価曲率．
	double Ry,			// in:	[m]:	径方向等価曲率．
	double bl_D,		// in:	[m]:	玉径
	double a,			// in:	[m]:	接触楕円長径．
	double b,			// in:	[m]:	接触楕円短径．
	double w,			// in:	[N]:	接触荷重．
	double E,			// in:	[Pa]:	物体0,1の等価ヤング率．
	double um,			// in:	[m/s]:	平均速度．
	double fratio,		// in:	[-]:	油膜荷重割合．
	double eta,			// in:	[Pas]:	油の粘度．
	double alpha,		// in:	[1/Pa]:	油の圧力粘度係数．
	double beta,		// in:	[1/K]:	油の温度粘度係数．
	double k,			// in:	[W/m2K]:油の熱伝導率．
	double xio			// in:	[m]:	メニスカス長さ
) {
	double p     = 1.5 * w / (Numeric::pi * a * b);
	double kap   = b / a;
	double beta_ = b * E / (Numeric::pi * p * Rx);
	double W     = 2 * Numeric::pi * p * b * b / (3 * kap * E * Rx * Rx);
	double U     = 2 * eta * um / E / Rx;
	double M     = pow(beta_, 0.25) * kap * W * pow(U, -0.75);
	double T     = 2 * a / Rx * ((0.77 * pow(beta_, -1/3) * pow(kap, 0.12) * pow(U, -1./12.) -1.8 * pow(beta_, -0.25)) / (1 + M / 6.6) + 1.8 * pow(beta_, -0.25)) * pow(U, 0.75);

	double t     = T * E * Rx * Rx * Rx;
	return t;
}

// EHL rolling resistance formula in Aihara's paper [ASME J. Trib., 109 (1987), 471] の式
double Tribology::AiharaR::calc
(						// out:	[Nm]:	トルク．
	double Rx,			// in:	[m]:	転がり方向等価曲率．
	double Ry,			// in:	[m]:	径方向等価曲率．
	double bl_D,		// in:	[m]:	玉径
	double a,			// in:	[m]:	接触楕円長径．
	double b,			// in:	[m]:	接触楕円短径．
	double w,			// in:	[N]:	接触荷重．
	double E,			// in:	[Pa]:	物体0,1の等価ヤング率．
	double um,			// in:	[m/s]:	平均速度．
	double fratio,		// in:	[-]:	油膜荷重割合．
	double eta,			// in:	[Pas]:	油の粘度．
	double alpha,		// in:	[1/Pa]:	油の圧力粘度係数．
	double beta,		// in:	[1/K]:	油の温度粘度係数．
	double k,			// in:	[W/m2K]:油の熱伝導率．
	double xio			// in:	[m]:	メニスカス長さ
) {
	double U = eta * um / E / Rx;
	double G = alpha * E;
	double W = w / E / Rx / Rx;
	//double W = 0.75 * w / E / Rx / a;										// 線接触の荷重パラメータの求め方に合わせる
	//double lc = a;// Tribology::RollerLength(a, Rx, Ry, xio);				// ころの有効長さ
	double L = eta * beta * um * um / k;
	double M = 176 / (1 + 0.29 * pow(L, 0.78)) / alpha * pow(G * U, 0.658) * pow(W, 0.31) * Rx * Rx * a;
	//double M = 176 / (1 + 0.29 * pow(L, 0.78)) / alpha * pow(G * U, 0.658) * pow(W, 0.31) * Rx * Rx * lc;
	return M;
}

// NTN TECHNICAL REVIEW No.82（2014）より，藤原の式．
double Tribology::Fujiwara::calc
(						// out:	[Nm]:	トルク．
	double Rx,			// in:	[m]:	転がり方向等価曲率．
	double Ry,			// in:	[m]:	径方向等価曲率．
	double bl_D,		// in:	[m]:	玉径
	double a,			// in:	[m]:	接触楕円長径．
	double b,			// in:	[m]:	接触楕円短径．
	double w,			// in:	[N]:	接触荷重．
	double E,			// in:	[Pa]:	物体0,1の等価ヤング率．
	double um,			// in:	[m/s]:	平均速度．
	double fratio,		// in:	[-]:	油膜荷重割合．
	double eta,			// in:	[Pas]:	油の粘度．
	double alpha,		// in:	[1/Pa]:	油の圧力粘度係数．
	double beta,		// in:	[1/K]:	油の温度粘度係数．
	double k,			// in:	[W/m2K]:油の熱伝導率．
	double xio			// in:	[m]:	メニスカス長さ
) {
	double U   = eta * um / E / Rx;
	double G   = alpha * E;
	double W   = w / E / Rx / Rx;
	double kap = b / a;
	double FHD = 44.6 * pow(U, 0.694) * pow(G, 0.961) * pow(W, 0.401) * (1 - 0.962 * exp(-0.818 * kap));
	double Fr  = FHD * Rx * Rx / alpha;

	double M  = Fr * Rx;
	return M;
}

// Baltacに実装されている"相原の式"．論文掲載の相原の式とと同様にGoksemの式を改良して作成したものだが，玉軸受用にチューニングされている（研究報告第457号に記載）
double Tribology::GoksemAihara::calc
(						// out:	[Nm]:	トルク．
	double Rx,			// in:	[m]:	転がり方向等価曲率．
	double Ry,			// in:	[m]:	径方向等価曲率．
	double bl_D,		// in:	[m]:	玉径
	double a,			// in:	[m]:	接触楕円長径．
	double b,			// in:	[m]:	接触楕円短径．
	double w,			// in:	[N]:	接触荷重．
	double E,			// in:	[Pa]:	物体0,1の等価ヤング率．
	double um,			// in:	[m/s]:	平均速度．
	double fratio,		// in:	[-]:	油膜荷重割合．
	double eta,			// in:	[Pas]:	油の粘度．
	double alpha,		// in:	[1/Pa]:	油の圧力粘度係数．
	double beta,		// in:	[1/K]:	油の温度粘度係数．
	double k,			// in:	[W/m2K]:油の熱伝導率．
	double xio			// in:	[m]:	メニスカス長さ
) {

	double lc =  Tribology::RollerLength(a, Rx, Ry, xio);		// ころの有効長さ（の半分?）
					
	double w_mj = 0.75 * w / a;			// 線荷重[N/m]
	double Wehl = w_mj / (E * Rx);		// 荷重パラメータ[-]
	double G = alpha * E;
	double U = eta * um / E / Rx;

	double L_GH = eta * beta * um * um / k;			// 熱負荷係数
	double D_GH = 11.81220746 * G * U / pow(Wehl, 1.5);			// 動負荷係数
	double Kt = exp((0.658228 + 0.0106383 * pow(L_GH, 0.3174603)) * log(D_GH) - (0.162625 + log(1 + 0.411483 * pow(L_GH, 0.707071))));
	
	double el = 2 * lc * 1000;					// ころの有効長さ[mm]
	double Da = bl_D * 1000;				// 玉径[mm]
	double wle = w / 9.8065 / el;
	double f_Q = pow((2 * fratio * wle / Da), 0.3);	// 荷重補正項（無次元数になるはずだが，なぜか無次元ではない）

	
	double Cst = this->Starvation(xio, D_GH, L_GH, b, Kt);	// スタベーション係数

	double tr = Cst * (1 / alpha) * fratio * Wehl * Kt * Rx;
	double Mr = tr * 2 * Rx * 2 * lc * f_Q;

	return Mr;
}

// ヒステリシス損失を計算．計算式は "角田和雄, 玉軸受の摩擦モーメントに関する研究（スラスト軸受の場合）機械学会誌27巻178号"に掲載
double Tribology::Kakuta::calc
(					// out: [Nm]:	ヒステリシス損失
	double fh,		// in:	[-]:	ヒステリシス係数
	double b,		// in:	[m]:	接触楕円短半径
	double Q		// in:	[N]:	転動体荷重
) {
	double Mrh = 0.375 * fh * b* Q;
	return Mrh;
}

// ヒステリシス損失がないとき
double Tribology::HysteresisNothing::calc
(					// out: [Nm]:	ヒステリシス損失
	double fh,		// in:	[-]:	ヒステリシス係数
	double b,		// in:	[m]:	接触楕円短半径
	double Q		// in:	[N]:	転動体荷重
) {
	return 0;
}

// Goksemの式のスタベーション係数
double Tribology::GoksemAihara::Starvation(
	double xio,		//	in:	[m]:	メニスカス長さ
	double D_GH,	//	in: [-]:	動負荷係数
	double L_GH,	//	in: [-]:	熱負荷係数
	double b,		//	in:	[m]:	接触楕円短半径
	double Kt		//	in:	[-]:	??
){
	// 十分潤滑 or スタベーション係数直接入力のときは補正無し
	if (xio <= 0) {
		return 1;
	}
	double Ca = 950 * pow(D_GH, 0.05) / (1000 + (0.2 * D_GH) * pow(L_GH, (0.66 - 0.003 * D_GH))); // スタベーションの影響[-]
	double Cb = pow(5 + D_GH, 0.068) * (pow(1 + L_GH, 0.06 / pow(D_GH, 0.3) - 0.035) - 0.05 * pow(L_GH, 1e-4));
	double mst = 1 + xio / b;
	double Kai = mst * sqrt(mst * mst - 1) - log(mst + sqrt(mst * mst - 1));
	double fai = pow(3 * Kai / (42 * Kt), 0.3);
	double Cst = Ca * pow(fai, Cb) / (1 + Ca * pow(fai, Cb));

	return Cst;

}


// 転がり摩擦がない状態を定義．
double Tribology::RollingResistanceNothing::calc(double Rx, double Ry, double bl_D, double a, double b, double w, double E, double um, double fratio, double eta, double alpha, double beta, double k, double xio) {
	return 0.0;
}



// HamrockDowsonの中央油膜厚さの式．発熱による粘度変化を考慮したい場合は下のErtelGrubin式をかけてください．（線形計算なのでそのまま掛け算できます．）（出典：Hamrock, Bernard J., and D. Dowson. "Isothermal elastohydrodynamic lubrication of point contacts. IV-Starvation results." (1976).）
double Tribology::HamrockDowsonHc::calc
(						// out:	[m]:	油膜厚さ．
	double w,			// in:	[N]:	接触荷重．
	double E,			// in:	[Pa]:	等価ヤング率．
	double Rx,			// in:	[m]:	転がり方向曲率半径．
	double Ry,			// in:	[m]:	周方向曲率半径．
	double alpha,		// in:	[1/Pa]:	圧力粘度係数．
	double eta,			// in:	[Pas]:	油の粘度．
	double u,			// in:	[m/s]:	転がり速度．
	double lm,			// in:	[m]:	メニスカス長さ．（<-1:自動計算，1~0:そのまま掛け算，>0:メニスカス長さ使用．）
	double bh			// in:	[m]:	接触楕円短径．
) {
	double W = w / E / Rx / Rx;
	double G = alpha * E;
	double U = eta * u / E / Rx;
	double k_ = 1.0339 * pow(Ry / Rx, 0.636);

	double H = 2.69 * pow(U, 0.67) * pow(G, 0.53) * pow(W, -0.067) * (1.0 - 0.61 * exp(-0.73 * k_));
	double h =  H * Rx;

	h *=  this->Starvation(H, Rx, lm, bh);
	return h;
}

// HamrockDowsonの最小油膜厚さの式．発熱による粘度変化を考慮したい場合は下のErtelGrubin式をかけてください．（線形計算なのでそのまま掛け算できます．）（出典：Hamrock, Bernard J., and D. Dowson. "Isothermal elastohydrodynamic lubrication of point contacts. IV-Starvation results." (1976).）
double Tribology::HamrockDowsonHmin::calc
(						// out:	[m]:	油膜厚さ．
	double w,			// in:	[N]:	接触荷重．
	double E,			// in:	[Pa]:	等価ヤング率．
	double Rx,			// in:	[m]:	転がり方向曲率半径．
	double Ry,			// in:	[m]:	周方向曲率半径．
	double alpha,		// in:	[1/Pa]:	圧力粘度係数．
	double eta,			// in:	[Pas]:	油の粘度．
	double u,			// in:	[m/s]:	転がり速度．
	double lm,			// in:	[m]:	メニスカス長さ．（<-1:自動計算，1~0:そのまま掛け算，>0:メニスカス長さ使用．）
	double bh			// in:	[m]:	接触楕円短径．
) {
	double W = w / E / Rx / Rx;
	double G = alpha * E;
	double U = eta * u / E / Rx;
	double k_ = 1.0339 * pow(Ry / Rx, 0.636);

	double H = 3.63 * pow(U, 0.68) * pow(G, 0.49) * pow(W, -0.073) * (1.0 - exp(-0.68 * k_));
	double h = H * Rx;

	h *= this->Starvation(H, Rx, lm, bh);
	return h;
}



// HamrockDowson中央油膜厚さのスタベーション係数．
double Tribology::HamrockDowsonHc::Starvation
(						// out:	[-]:	スタベーション係数．
	double H,			// in:	[-]:	油膜パラメタ．
	double Rx,			// in:	[m]:	転がり方向曲率半径．
	double lm,			// in:	[m]:	メニスカス長さ．（<-1:自動計算，1~0:そのまま掛け算，>0:メニスカス長さ使用．）
	double bh			// in:	[m]:	接触楕円短径．
) {
	if (lm <=  -1.0)
		return 1.0;

	if (lm <=  0.0)
		return fabs(lm);

	double f  = lm / bh;
	double fc = 3.06 * pow(H * Rx * Rx / bh / bh, 0.58);
	if (f > fc)
		return 1.0;

	return pow(f / fc, 0.29);
}

// HamrockDowson最小油膜厚さのスタベーション係数．
double Tribology::HamrockDowsonHmin::Starvation
(						// out:	[-]:	スタベーション係数．
	double H,			// in:	[-]:	油膜パラメタ．
	double Rx,			// in:	[m]:	転がり方向曲率半径．
	double lm,			// in:	[m]:	メニスカス長さ．（<-1:自動計算，1~0:そのまま掛け算，>0:メニスカス長さ使用．）
	double bh			// in:	[m]:	接触楕円短径．
) {
	if (lm <= -1.0)
		return 1.0;

	if (lm <= 0.0)
		return fabs(lm);

	double f = lm / bh;
	double fc = 3.34 * pow(H * Rx * Rx / bh / bh, 0.58);
	if (f > fc)
		return 1.0;

	return pow(f / fc, 0.25);
}

// 油膜厚さがない状態を定義．
double Tribology::FilmThicknessNothing::calc(double w, double E, double Rx, double Ry, double alpha, double eta, double u, double lm, double bh) {
	return 0.0;
}

// 発熱による粘度変化を考慮した補正係数．
double Tribology::ErtelGrubin
(						// out:	[-]:	補正係数（0~1の実数）．油膜厚さにかけ算を行う．
	double eta,			// in:	[Pas]:	油の粘度．
	double beta,		// in:	[1/K]:	温度粘度係数．
	double k,			// in:	[W/m2K]:油の熱伝導率．
	double u			// in:	[m/s]:	転がり速度．
) {
	double L = eta * beta * u * u / k;
	double phi =  3.94 / (3.94 + pow(L, 0.62));
	return phi;
}

// 相原の式からトラクション係数を計算する．
double Tribology::AiharaT::calc
(						// out:	[-]:	トラクション係数．
	double eta0,		// in:	[Pas]:	油の粘度．
	double p0,			// in:	[Pa]:	最大面圧．
	double v0,			// in:	[m/s]:	滑り速度．
	double u0			// in:	[m/s]:	転がり速度．
) {
	if (eta0 <= 0.0)
		return 0.0;

	if (p0 <= 0.0)
		return 0.0;

	if (v0 <= 0.0)
		return 0.0;

	if (u0 <= 0.0)
		return 0.0;

	double s    = v0 / u0;
	double eta  = Unit::Pas2cP(eta0);
	double u    = Unit::ms2mms(u0);
	double pmax = Unit::Pa2kgfmm2(p0);
	double smax = 1.2448e7 * pow(eta, -1.5) * pow(u, -0.75) * pow(pmax, -1.5);

	double mumax;
	if (pmax > 94.404)
		mumax = 4.078e-2 * pow(eta, 0.1) * pow(u, -0.12) * pow(pmax, 0.18);
	else
		mumax = 1.008e-4 * pow(eta, 0.1) * pow(u, -0.12) * pow(pmax, 1.5);

	double mu = mumax * (1.0 - exp(-2.6 * s / smax));

	return mu;
}


// 部材の材料と接触部剛性から，減衰を考慮した接触力を求める．"15_非線形バネと衝突のモデル化"参考．
// 計算式は"日本惑星科学会誌, 2004, 13.4: 233-240."の非線形ばねへの手計算による適応． 
double Tribology::Tsuji::calc
(
	double k,		// in:	[N/x^1.5]:	非線形剛性．
	double zeta,	// in:	[-]		 :	減衰比．
	double m,		// in:	[kg]	 :	等価質量．
	double v,		// in:	[m/s]	 :	接近速度．
	double x		// in:	[m]		 :	接近量．
) {
	double c = 2 * zeta * sqrt(1.5 * m * k);
	double Fn = k * pow(x, 1.5) + c * v * pow(x, 0.25);
	// 負値のときは0に補正．
	return std::max(Fn, 1e-20);
}

// Kelvin & Voigtの接触モデルを用いて減衰を考慮した接触力を計算 
double Tribology::KelvinVoigt::calc
(
	double k,		// in:	[N/x]:	線形剛性．
	double zeta,	// in:	[-]		 :	減衰比．
	double m,		// in:	[kg]	 :	等価質量．
	double v,		// in:	[m/s]	 :	接近速度．
	double x		// in:	[m]		 :	接近量．
) {
	double c = 2 * zeta * sqrt(m * k);
	double Fn = k * x + c * v;
	// 負値のときは0に補正．
	return std::max(Fn, 1e-20);
}


// 旧動解析BRAINで実装されていたHamrockDowsonの油膜厚さの式．
//double Tribology::HamrockDowson_old::calc
//(						// out:	[m]:	油膜厚さ．
//	double w,			// in:	[N]:	接触荷重．
//	double E,			// in:	[Pa]:	等価ヤング率．
//	double Rx,			// in:	[m]:	転がり方向曲率半径．
//	double Ry,			// in:	[m]:	周方向曲率半径．
//	double alpha,		// in:	[1/Pa]:	圧力粘度係数．
//	double eta,			// in:	[Pas]:	油の粘度．
//	double u,			// in:	[m/s]:	転がり速度．
//	double lm,			// in:	[m]:	メニスカス長さ．（<-1:自動計算，1~0:そのまま掛け算，>0:メニスカス長さ使用．）
//	double bh			// in:	[m]:	接触楕円短径．
//) {
//	double W = w / E / Rx / Rx;
//	double G = alpha * E;
//	double U = eta * u / E / Rx;
//
//	double H = 3.63 * pow(U, 0.68) * pow(G, 0.49)
//		* pow(W, -0.073) * (1.0 - exp(-0.7031 * pow(Ry/Rx, 0.636)));
//
//	// Starvation
//	double starv;
//	if (lm <=  -1.0)
//		starv = 1.0;
//
//	else if (lm <=  0.0)
//		starv = fabs(lm);
//
//	else {
//		double f  = 1.0 + lm / bh;
//		double fc = 1.0 + 3.34 * pow(H * Rx * Rx / bh / bh, 0.56);
//		if (f > fc)
//			starv = 1.0;
//		else
//			starv = pow((f - 1.0)/(fc - 1.0), 0.25);
//	}
//	H *= starv;
//
//	double h  =  H * Rx;
//	return h;
//}
//
//// 旧動解析BRAINに実装されていた式．おそらく間違いだと思うが一応掲載．
//double Tribology::Aihara_old::calc
//(						// out:	[Nm]:	トルク．
//	double Rx,			// in:	[m]:	周方向等価曲率．
//	double Ry,			// in:	[m]:	転がり方向等価曲率．
//	double a,			// in:	[m]:	接触楕円長径．
//	double b,			// in:	[m]:	接触楕円短径．
//	double w,			// in:	[N]:	接触荷重．
//	double E,			// in:	[Pa]:	物体0,1の等価ヤング率．
//	double um,			// in:	[m/s]:	平均速度．
//	double eta,			// in:	[Pas]:	油の粘度．
//	double alpha,		// in:	[1/Pa]:	油の圧力粘度係数．
//	double beta,		// in:	[1/K]:	油の温度粘度係数．
//	double k			// in:	[W/m2K]:油の熱伝導率．
//) {
//	double W = w / E / Ry / Ry;
//	double G = alpha * E;
//	double U = eta * um / E / Ry;
//
//	double L = eta * beta * um * um / k;
//	double p = 3.94 / (3.94 + pow(L, 0.62));
//	double H = (3.63 * pow(U, 0.68) * pow(G, 0.49) * pow(W, -0.073)
//		* (1.0 - exp(-0.7031 * pow(Rx / Ry, 0.636)))) * p;
//
//	double D  = 11.81221 * G * U / W / sqrt(W);
//	double c  = exp((0.658228 + 0.0106383 * pow(L, 0.3174603)) * log(D)
//		- 0.162625 - log(1.0 + 0.411483 * pow(L, 0.707071)));
//	double tr = c * W * Ry / alpha;
//
//	double M = 2.0 * pow(w, 0.3) * tr * Ry * a;
//	return M;
//	//double tr = c * fratio * W * ry / alpha_;
//	// w *= fratio;
//	//double lm			// Meniscus length
//	//double starv;
//	//if (lm <= 0.0)
//	//	starv = 1.0;
//	//else {
//	//	double a1;
//	//	double f = 1.0 + lm / b;
//	//	if ((0.66-0.003*D) > 0.0)
//	//		a1 = 950.0 * pow(D, 0.05) / (1000.0+(0.2*D+9.0)*pow(L, 0.66-0.003*D));
//	//	else
//	//		a1 = 950.0 * pow(D, 0.05) / 1000.0;
//	//	double b  = (pow(1.0+L, 0.06*pow(D, -0.3)-0.035)-0.05*pow(L, 0.0001))*pow(5.0+D, 0.068);
//	//	double f2 = sqrt(f * f - 1.0);
//	//	double s  = (f2 * f - log(f + f2)) / c;
//	//	if (s>0.0) {
//	//		s = 0.6551853*pow(s, 2.0/3.0);
//	//		double s2 = a1*pow(s, b);
//	//		starv = s2 / (1.0 + s2);
//	//	}
//	//	else {
//	//		starv = 0.0;
//	//	}
//	//}
//	//double M = starv * 2.0 * pow(w, 0.3) * tr * ry * a;
//	//if (fratio <= 0.0)
//	//	M = 0.0;
//};
//

