#include "pch.h"


class TribologyTest : public ::testing::Test {
};

// スライス荷重の総和が元の荷重と一致を確認．
TEST_F(TribologyTest, ForceSlice0) {

	double fs[100];
	Tribology::SliceForceRatio(100, fs);
	double return_ = 0.0;
	for (int i = 0; i < 100; i++)
		return_ += fs[i];

	double answer_ = 1.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);

	return;
};

// 4分割時のスライス荷重を計算．定積分を解析的に解いた手計算と比較．
TEST_F(TribologyTest, ForceSlice1) {

	double fs[4];
	Tribology::SliceForceRatio(4, fs);
	double return_0 = fs[0] * 100;
	double answer_0 = 15.625;	// = 5/32 
	double delta_0 = abs(return_0 - answer_0);
	EXPECT_GT(1e-2, delta_0);

	double return_1 = fs[1] * 100;
	double answer_1 = 34.375;	// = 11/32
	double delta_1 = abs(return_1 - answer_1);
	EXPECT_GT(1e-2, delta_1);

	return;
};

// スライス荷重の総和が元の荷重と一致を確認．
TEST_F(TribologyTest, SliceForceRatio2d0) {

	double fs[10];
	Tribology::SliceForceRatio2d(10, 4, fs);
	double return_ = 0.0;
	for (int i = 0; i < 10; i++)
		return_ += fs[i];

	double answer_ = 0.25;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);

	return;
};

// 等価ヤング率のテスト．手計算との一致を確認．
TEST_F(TribologyTest, ReducedYoung0) {

	double return_ = Tribology::ReducedYoung(10.0, 0.3, 10.0, 0.3);
	double answer_ = 10.99;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);

	return;
};

// 合成粗さのテスト．手計算との一致を確認．
TEST_F(TribologyTest, CompositeRoughness0) {

	double return_ = Tribology::CompositeRoughness(10.0, 10.0);
	double answer_ = 14.1421356;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);

	return;
};

// 換算質量のテスト．手計算との一致を確認．
TEST_F(TribologyTest, ReducedMass0) {

	double return_ = Tribology::ReducedMass(10.0, 10.0);
	double answer_ = 5.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);

	return;
};


// ラムダ値と推定される荷重支持割合の比較テスト．ラムダ値３で大体１になることを確認．
TEST_F(TribologyTest, ForceRatio0) {

	double ret0 = Tribology::ForceRatio(-1e-9);
	double ret1 = Tribology::ForceRatio(1e-9);
	double ret2 = Tribology::ForceRatio(3.0);

	EXPECT_GT(1e-2, abs(0 - ret0));
	EXPECT_GT(1e-2, abs(0 - ret1));
	EXPECT_GT(1e-2, abs(1 - ret2));

	return;
}

// 滑り率微小と極大で，それぞれ狙い通りの値となっているかテスト．
TEST_F(TribologyTest, Tangent0) {

	Tribology::Tangent TG;

	double mu = 0.10;
	double sigma0 = 1e-4 / 1e0;
	double mu0 = TG.calc(mu, 1e-4, 1e0, 1);
	EXPECT_GT(1e-2, abs(mu0 - 0.6366 * sigma0 * mu) / mu0);
	
	double sigma1 = 1e-1 / 1e0;
	double mu1 = TG.calc(mu, 1e-1, 1e0, 1e3);
	EXPECT_GT(1e-2, abs(mu1 - mu) / mu1);

	return;
}

// 運動力学IIにて値の妥当性を確認．
TEST_F(TribologyTest, Hertz0) {

	double Rx = 1 / 0.180868e3;
	double Ry = 1 / 0.005273e3;
	double E  = 228.57142857142856e9;
	double dx = 9.164e-6;

	Tribology::BrewHamrock HZ;
	double F, k, a, b;
	HZ.calc(Rx, Ry, dx, E, k, a, b);
	F = k * pow(dx, 1.5);
	double F_ = 972.3894714547508;
	double k_ = 35051985801.205666;
	double a_ = 0.0016211389380570802;
	double b_ = 0.00016553444618622523;

	EXPECT_GT(1e-2, abs(F - F_) / abs(F_));
	EXPECT_GT(1e-2, abs(k - k_) / abs(k_));
	EXPECT_GT(1e-2, abs(a - a_) / abs(a_));
	EXPECT_GT(1e-2, abs(b - b_) / abs(b_));

	return;
};

// 図解トライボロジー摩擦の科学と潤滑 p.32．にて値の妥当性を確認．
TEST_F(TribologyTest, Hertz1) {

	double Rx = 10.7e-3;
	double Ry = 20e-3;
	double E  = 231e9;
	double dx = 1.89e-6;

	Tribology::BrewHamrock HZ;
	double F, k, a, b;
	HZ.calc(Rx, Ry, dx, E, k, a, b);
	F = k * pow(dx, 1.5);
	double F_ = 50;
	double a_ = 0.208e-3;
	double b_ = 0.135e-3;

	EXPECT_GT(1e-2, abs(F - F_) / abs(F_));
	EXPECT_GT(1e-2, abs(a - a_) / abs(a_));
	EXPECT_GT(1e-2, abs(b - b_) / abs(b_));

	return;
};

// 転がり摩擦の異なる3社から出されたそれぞれの式．藤原の式が論文と一致していること，また3つの式でオーダーが一致していることを確認．
TEST_F(TribologyTest, RollingResistance0) {

	double D  = 4.3650e-3;
	double R  = 40.635e-3;
	double r  = 2.3571e-3;
	double E  = 207.00e+9;
	double nu = 0.3000e+0;
	double W  = 136.00e+0;
	double wIR= 10000 * Numeric::pi / 30;
	double dx = 4.5674e-6;
	double a  = 429.50e-6;
	double b  = 76.762e-6;
	double eta= 43.691E-3;
	double Rx = 1 / (2 / D + 1 / R);
	double Ry = 1 / (2 / D - 1 / r);
	double Eq = 2.0 / ((1.0 - nu * nu) / E + (1.0 - nu * nu) / E);
	double Dpcd = 2 * R + D;
	double wBL= (Dpcd / D - D / Dpcd) * 0.5 * wIR;
	double u  = D / 2 * wBL;
	double alpha = 0.2 / 9.8067e6;
	double beta  = 7.5e-4;
	double k     = 0.145;
	double fratio = 0, lm = -1;// 使わない
	Tribology::AiharaR AH;
	double T_AH = AH.calc(Rx, Ry, D, a, b, W, Eq, u, fratio, eta, alpha, beta, k, lm);

	Tribology::Houpert HP;
	double T_HP = HP.calc(Rx, Ry, D, a, b, W, Eq, u, fratio, eta, alpha, beta, k, lm);

	Tribology::Fujiwara FW;
	double T_FW = FW.calc(Rx, Ry, D, a, b, W, Eq, u, fratio, eta, alpha, beta, k, lm);

	double T_AH_ = 0.4815508324699387e-3;
	double T_HP_ = 0.7780845243544707e-3;
	double T_FW_ = 0.2931298106527921e-3;

	EXPECT_GT(1e-2, abs(T_AH - T_AH_));
	EXPECT_GT(1e-2, abs(T_HP - T_HP_));
	EXPECT_GT(1e-2, abs(T_FW - T_FW_));

	return;
};

// baltacに実装されている転がり粘性抵抗の式と計算結果が同じになっているか確認．
TEST_F(TribologyTest, RollingResistance1) {
	double D = 11.1125000000000e-3; // 玉径[m]
	double Rx = 6.34651992397817 * 1e-3;	// 等価半径[m](転がり方向)
	double a = 0.320081205545215 * 1e-3;	// 長半径[m]
	double b = 5.418495251937330e-5;	// 短半径[m]
	double W = 1.84083129679812 * 9.8065;	// 荷重[N]
	double Eq = 23157.5499508680 * 9.8065 * 1e6;	// 等価ヤング率[Pa]

	double Ry = 98.0563663691003e-3;	// 等価半径[m](軸方向)
	double u = 17.3140034332462e-3;		// 転がり速度[m/s]
	double alpha = 0.21 / 9.8067e6;		// 圧力粘度係数
	double eta = 7.579516957245301e-9 * 9.8065 * 1e6;	// 粘性[kgf/mm^2・s]⇒[Pa・s]
	double dvdt = 4.309343950426442e-10;// ???
	double beta = 0.056855126;// = dvdt / eta;	// 温度粘度係数[1/K]
	double k = 1.477550000000000e-002 * 9.8065;	// 熱伝導率[kgf/s・K]⇒[W/(m・K)]
	double lm = -1;
	double Tr_baltac = 9.029644084675681626e-6;	// 転がり粘性抵抗[kgfmm](baltac計算値)
	double fratio = 0.433752349655000;
	Tribology::GoksemAihara GA;
	double Tr_norm = GA.calc(Rx, Ry, D, a, b, W, Eq, u, fratio, eta, alpha, beta, k, lm);
	EXPECT_NEAR(Tr_norm, Tr_baltac, Tr_baltac * 1e-2);
	// デバッグで以下の数値になっていることを確認することが望ましい
	// D_GH = 0.322606730797286;
	// L_GH = 8.743070312928715e-006;
	// w_mj = 42,298.903682820381025;
	// GG = 4863.08548968228;
	// U = 8.929170740443306e-013
	// tr = 1.52212416932625653725;[N/m]
	// W = 2.934858572384527e-005;
	// Kt = 0.403451221364356
	// f_Q = 0.604784692973623
	// Mr = 9.029644084675681626e-6
	return;
};

// 油膜厚さの式．手計算（\docs\00_Theory\08. 点接触の油膜厚さ）との一致を確認．
TEST_F(TribologyTest, HamrockDowson_case1) {

	double D  = 4.3650e-3;
	double R  = 40.635e-3;
	double r  = 2.3571e-3;
	double E  = 207.00e+9;
	double nu = 0.3000e+0;
	double W  = 136.00e+0;
	double wIR= 10000 * Numeric::pi / 30;
	double a  = 429.50e-6;
	double b  = 76.762e-6;
	double eta= 43.691E-3;
	double Rx = 1 / (2 / D + 1 / R);
	double Ry = 1 / (2 / D - 1 / r);
	double Eq = 2.0 / ((1.0 - nu * nu) / E + (1.0 - nu * nu) / E);
	double Dpcd = 2 * R + D;
	double wBL= (Dpcd / D - D / Dpcd) * 0.5 * wIR;
	double u  = D / 2 * wBL;
	double alpha = 0.2 / 9.8067e6;

	Tribology::HamrockDowsonHc HD;

	// (a) 枯渇潤滑の時
	double h0 = HD.calc(W, Eq, Rx, Ry, alpha, eta, u, b, b);
	double h0_ = 1.09573218705714E-06;

	EXPECT_NEAR(h0, h0_, 1e-8);

	// (b) 十分潤滑の時
	double h1 = HD.calc(W, Eq, Rx, Ry, alpha, eta, u, 3 * b, b);
	double h1_ = 1.33431846824466E-06;

	EXPECT_NEAR(h1, h1_, 1e-8);

	return;
};

// 油膜厚さの式．手計算（\docs\00_Theory\08. 点接触の油膜厚さ）との一致を確認．
TEST_F(TribologyTest, HamrockDowson_case2) {

	double D  = 4.3650e-3;
	double R  = 40.635e-3;
	double r  = 2.3571e-3;
	double W  = 100.00e+0;
	double a  = 400.00e-6;
	double lm = -1.0;
	double b  = 80.00e-6;
	double eta= 40.0E-3;
	double Rx = 20.0e-3;
	double Ry = 30.0e-3;
	double Eq = 200.00e+9;
	double u = 10;
	double alpha = 2.0e-8;
	 
	Tribology::HamrockDowsonHmin HD;

	// (a) 枯渇潤滑の時
	double h0 = HD.calc(W, Eq, Rx, Ry, alpha, eta, u, lm, b);
	double h0_ = 1.07933E-06;
	EXPECT_NEAR(h0, h0_, 1e-8);
	return;
};

// HamrockDowsonのスタベーション係数の確認．
TEST_F(TribologyTest, Starvation_case1) {

	double H  = 0.00064420825708110159;
	double Rx = 0.0020712532842879663;
	double bh = 7.6761999999999995e-05;
	Tribology::HamrockDowsonHc HD;
	
	// (a) メニスカス長さが -1 以下であれば，1を出力．
	double s0 = HD.Starvation(H, Rx, -1.0, bh);
	EXPECT_NEAR(s0, 1.0, 1e-3);

	// (b) メニスカス長さが -1<lm<0 であれば，その絶対値を補正係数とする
	double s1 = HD.Starvation(H, Rx, -0.3, bh);
	EXPECT_NEAR(s1, 0.3, 1e-3);

	// (c) メニスカス長さが 0を超える場合，計算値を用いる．
	double s2 = HD.Starvation(H, Rx, bh, bh);
	EXPECT_NEAR(s2, 0.821192401315268, 1e-3);

	// (d) メニスカス長さが 0を超える場合，計算値を用いる．
	//     ただし計算値が1を超える場合，1を出力
	double s3 = HD.Starvation(H, Rx, bh * 3, bh);
	EXPECT_NEAR(s3, 1, 1e-3);

	return;
};

// HamrockDowsonのスタベーション係数の確認．
TEST_F(TribologyTest, Starvation_case2) {

	double H = 6e-4;
	double Rx = 2e-3;
	double bh = 6e-5;
	Tribology::HamrockDowsonHmin HD;

	// (a) メニスカス長さが -1 以下であれば，1を出力．
	double s0 = HD.Starvation(H, Rx, -1.0, bh);
	EXPECT_NEAR(s0, 1.0, 1e-3);

	// (b) メニスカス長さが -1<lm<0 であれば，その絶対値を補正係数とする
	double s1 = HD.Starvation(H, Rx, -0.3, bh);
	EXPECT_NEAR(s1, 0.3, 1e-3);

	// (c) メニスカス長さが 0を超える場合，計算値を用いる．
	double s2 = HD.Starvation(H, Rx, bh, bh);
	EXPECT_NEAR(s2, 0.784506626, 1e-3);

	// (d) メニスカス長さが 0を超える場合，計算値を用いる．
	//     ただし計算値が1を超える場合，1を出力
	double s3 = HD.Starvation(H, Rx, bh * 3, bh);
	EXPECT_NEAR(s3, 1, 1e-3);

	return;
};

// HamrockDowsonのスタベーション係数の確認．
TEST_F(TribologyTest, ErtelGrubin0) {

	double eta  = 43.691E-3;
	double beta = 7.5e-4;
	double k    = 0.145;
	double u    = 22.36094193;

	double phi = Tribology::ErtelGrubin(eta, beta, k, u);
	EXPECT_GT(1e-2, abs(phi- 0.938371825085712) / phi);

	return;
};

// 相原トラクション係数の式の条件分岐の確認．
TEST_F(TribologyTest, AiharaT0) {

	Tribology::AiharaT AT;
	EXPECT_GT(1e-2, AT.calc(-1.0, 1.0, 1.0, 1.0));
	EXPECT_GT(1e-2, AT.calc(1.0, -1.0, 1.0, 1.0));
	EXPECT_GT(1e-2, AT.calc(1.0, 1.0, -1.0, 1.0));
	EXPECT_GT(1e-2, AT.calc(1.0, 1.0, 1.0, -1.0));

	return;
};

// 相原の式からトラクション係数の確認(荷重大)
TEST_F(TribologyTest, AiharaT1) {

	Tribology::AiharaT AT;
	double eta0 = 72.9029224284153 * 0.001; // 油の粘度 [cP→Pas]
	double p0 = 149.980489558372 * 9.80665 * 1e6;// 最大面圧[kgf/mm^2→Pa]
	double v0 = 48.8863222878130 / 1000; // 滑り速度[mm/s → m/s]
	double u0 = 896.117562552404 / 1000; // 転がり速度[mm/s → m/s]
	double mu = AT.calc(eta0, p0, v0, u0);	// トラクション係数（実際）
	double mu_ex = 6.017253174237477e-002;	// トラクション係数（baltacデバッグの計算値）
	ASSERT_NEAR(mu_ex, mu, 1e-4);
	return;
};


// 相原の式からトラクション係数の確認(荷重小)
TEST_F(TribologyTest, AiharaT2) {

	Tribology::AiharaT AT;
	double eta0 = 72.9029224284153 * 0.001; // 油の粘度 [cP→Pas]
	double p0 = 70.8668300974355 * 9.80665 * 1e6;// 最大面圧[kgf/mm^2→Pa]
	double v0 = 18.7016744160132 / 1000; // 滑り速度[mm/s → m/s]
	double u0 = 904.002621182389 / 1000; // 転がり速度[mm/s → m/s]
	double mu = AT.calc(eta0, p0, v0, u0);	// トラクション係数（実際）
	double mu_ex = 9.483619302366691e-003;	// トラクション係数（baltacデバッグの計算値）
	ASSERT_NEAR(mu_ex, mu, 1e-4);
	return;
};


// 0を返すはずのメソッドがちゃんと0を返すかのテスト．
TEST_F(TribologyTest, Nothing0) {

	Tribology::CoulombNothing n0;
	EXPECT_GT(1e-2, n0.calc(1,1,1,1));
	
	Tribology::RollingResistanceNothing n1;
	EXPECT_GT(1e-2, n1.calc(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1));

	Tribology::FilmThicknessNothing n2;
	EXPECT_GT(1e-2, n2.calc(1, 1, 1, 1, 1, 1, 1, 1, 1));

	return;
}


// (1) 非線形減衰(Tsuji)を手計算と比較．
TEST_F(TribologyTest, Tsuji_test1) {
	double m = 1;
	double k = 1;
	double zeta = 0.2;
	double v = 1.0;
	double x = 16 * 1e-8;
	Tribology::DampingForce *DF = new Tribology::Tsuji();
	double Fc = DF->calc(k, zeta, m, v, x);
	double Fc_ex = 0.489897948*0.02;
	ASSERT_NEAR(Fc, Fc_ex, 1e-6);

	return;
}

// (2) 非線形減衰(Tsuji)で接触力が負の時0に補正するか確認．
TEST_F(TribologyTest, Tsuji_test2) {
	double m = 1;
	double k = 1;
	double zeta = 0.2;
	double v = -10.0;
	double x = 16 * 1e-8;
	Tribology::DampingForce *DF = new Tribology::Tsuji();
	double Fc = DF->calc(k, zeta, m, v, x);
	double Fc_ex = 0;
	ASSERT_NEAR(Fc, Fc_ex, 1e-6);

	return;
}

// (1) 線形減衰(KelvinVoigt)を手計算と比較．
TEST_F(TribologyTest, KelvinVoigt_test1) {
	double m = 1;
	double k = 1;
	double zeta = 0.2;
	double v = 1.0;
	double x = 0.5;
	Tribology::DampingForce *DF = new Tribology::KelvinVoigt();
	double Fc = DF->calc(k, zeta, m, v, x);
	double Fc_ex = 0.9;
	ASSERT_NEAR(Fc, Fc_ex, 1e-4);
	return;
}

// (2) 線形減衰(KelvinVoigt)で接触力が負の時0に補正するか確認．
TEST_F(TribologyTest, KelvinVoigt_test2) {
	double m = 1;
	double k = 1;
	double zeta = 0.2;
	double v = -2.0;
	double x = 0.5;
	Tribology::DampingForce *DF = new Tribology::KelvinVoigt();
	double Fc = DF->calc(k, zeta, m, v, x);
	double Fc_ex = 0;
	ASSERT_NEAR(Fc, Fc_ex, 1e-4);
	return;
}


// 接近量が微小量の時，減衰ありと無しで面圧がどのくらい変わるかを確認
TEST_F(TribologyTest, Pmax_test1) {
	double Rx = 5e-3;
	double dx = 1e-6;
	double E = 200e9;
	double zeta = 0.3;
	double m = Rx * Rx * Rx * Numeric::pi * 1.3333333333333333333333333333333;
	double v = 0.1;

	double k, a, b;
	Tribology::DampingForce *DF = new Tribology::KelvinVoigt();
	Tribology::Hertz *HZ = new Tribology::BrewHamrock();
	// BrewHamrock の近似式から剛性，接触楕円を求める． 
	HZ->calc(Rx, Rx, dx, E, k, a, b);
	double Fc = DF->calc(k, zeta, m, v, dx);
	std::cout << Fc << "\t" << k * pow(dx, 1.5) << std::endl;
	double Fc_ex = 0;
	ASSERT_NEAR(Fc, Fc_ex, 1e-4);
	return;
}