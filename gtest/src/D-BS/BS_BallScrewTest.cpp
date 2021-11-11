/*******************************************************************************
!						"BS_BallScrewTest.cpp"
!														2021/05/19	[Core-T] 柏原
!	初期値計算のメソッドをメインにテスト．入力・出力用の関数は正しく機能するものとして考える．
!	ただし，玉初期座標など若干複雑な計算をしているものはテストする．
!
!*******************************************************************************/

#include "pch.h"
#include "BS_FileIn_stab.h"
#include "BS_FileOut_stab.h"
using Eigen::IOFormat;
using namespace std;

#define cos30 0.86602540378443864676372317075294
#define sin30 0.5

class BS_BallSrewTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	BS_BallScrew BS;
	BS_FileIn_stab IN;
	BS_FileOut_stab OUT;

	virtual void SetUp() {
	}

	// ボールねじの幾何形状を定義
	// 溝一つ，螺旋の傾きが30°の単純形状にする
	void SingleNut_init() {
		// 潤滑理論
		IN.tribology.rollingresistance = BS_In::Tribology::RollingResistance::RollingResistanceNothing;
		IN.tribology.coulomb = BS_In::Tribology::Tangent;
		IN.tribology.filmThickness = BS_In::Tribology::FilmThicknessNothing;
		IN.tribology.coulomb_slope = 1000;

		// 回路・玉
		this->IN.circuit.resize(1);
		this->IN.circuit[0].th0 = 0;
		this->IN.circuit[0].th1 = Numeric::pi * 2;
		this->IN.circuit[0].is = 0;
		this->IN.circuit[0].inut = 0;
		this->IN.circuit[0].ball.resize(5);
		for (int i = 0; i < 5; i++) {
			this->IN.circuit[0].ball[i].density = 7830;
			this->IN.circuit[0].ball[i].young = 207760000000;
			this->IN.circuit[0].ball[i].poisson = 0.3;
			this->IN.circuit[0].ball[i].rms = 0.00000002;
			this->IN.circuit[0].ball[i].r = 10e-3;
		}

		// シャフト物性値(螺旋の初期位相角は0°)
		this->IN.shaft.resize(1);
		this->IN.shaft[0].m = 1;
		this->IN.shaft[0].young = 207760000000;
		this->IN.shaft[0].poisson = 0.3;
		this->IN.shaft[0].Ix = 1, IN.shaft[0].Iyz = 1, IN.shaft[0].Iyz = 1;

		this->IN.shaft[0].spiral.resize(1);
		this->IN.shaft[0].spiral[0].groove[0].eta[1] = 0;
		this->IN.shaft[0].spiral[0].groove[1].eta[1] = 0;
		this->IN.shaft[0].spiral[0].groove[0].r = 11e-3;
		this->IN.shaft[0].spiral[0].groove[1].r = 1;
		// 溝R中心PCD[m]（溝pcd半径，0.0357535 m くらい）
		this->IN.shaft[0].spiral[0].groove[0].eta[0] = -0.5e-3;
		this->IN.shaft[0].spiral[0].groove[1].eta[0] = -0.5e-3;

		this->IN.shaft[0].spiral[0].groove[0].sigma = 0.00000006;
		this->IN.shaft[0].spiral[0].groove[1].sigma = 0.00000006;
		this->IN.shaft[0].spiral[0].alp = 0;
		this->IN.shaft[0].spiral[0].r = 50e-3;
		this->IN.shaft[0].spiral[0].l = 181.2879e-3;// (= 2πr * tan30°)

		// ナット物性値
		this->IN.nut.resize(1);
		this->IN.nut[0].spiral.resize(1);
		this->IN.nut[0].spiral[0].groove[0].eta[1] = 0;
		this->IN.nut[0].spiral[0].groove[1].eta[1] = 0;
		this->IN.nut[0].spiral[0].groove[0].r = 11e-3;
		this->IN.nut[0].spiral[0].groove[1].r = 1e3;
		this->IN.nut[0].spiral[0].groove[0].eta[0] = 0.5e-3;
		this->IN.nut[0].spiral[0].groove[1].eta[0] = 0.5e-3;
		this->IN.nut[0].spiral[0].groove[0].sigma = 0.00000006;
		this->IN.nut[0].spiral[0].groove[1].sigma = 0.00000006;
		this->IN.nut[0].spiral[0].alp = 0;
		this->IN.nut[0].spiral[0].r = 50e-3;
		this->IN.nut[0].spiral[0].l = 181.2879e-3;// (= 2πr * tan30°)

		this->IN.nut[0].m = 1;
		this->IN.nut[0].young = 207760000000;
		this->IN.nut[0].poisson = 0.3;
		this->IN.nut[0].Ix = 1, this->IN.nut[0].Iyz = 1, this->IN.nut[0].Iyz = 1;

		// Pairクラス物性値
		int n_ball = this->IN.circuit[0].ball.size();
		this->IN.BallShaftPair.resize(n_ball);
		for (int i = 0; i < n_ball; i++) {
			this->IN.BallShaftPair[0].groove[0].zeta = 0.2;	// 減衰項
			this->IN.BallShaftPair[0].groove[0].mu = 0.1;		// 摩擦係数
			this->IN.BallShaftPair[0].groove[1].zeta = 0.2;	// 減衰項
			this->IN.BallShaftPair[0].groove[1].mu = 0.1;		// 摩擦係数
		}
		this->IN.BallNutPair.resize(n_ball);
		for (int i = 0; i < n_ball; i++) {
			this->IN.BallNutPair[0].groove[0].zeta = 0.2;	// 減衰項
			this->IN.BallNutPair[0].groove[0].mu = 0.1;		// 摩擦係数
			this->IN.BallNutPair[0].groove[1].zeta = 0.2;	// 減衰項
			this->IN.BallNutPair[0].groove[1].mu = 0.1;		// 摩擦係数
		}

		this->IN.oil.eta = 0.5;	// 粘度（Pa*s）(適当)
		this->IN.oil.beta = 0.5;	// 温度粘度係数（適当）
		this->IN.oil.k = 0.145;	// 油熱伝導率[W/(m・K)]
		this->IN.oil.alpha = 0.2;	// 圧力粘度係数[mm2/kgf]
		this->IN.oil.lm = -1;	// メニスカス長さ

		// 荷重（必要に応じて各テストケースの中で定義）
		this->IN.load.resize(1);
		this->IN.load[0].x[0] = 0;
		this->IN.load[0].x[1] = 0;
		this->IN.load[0].x[2] = 0;
		this->IN.load[0].F[0] = 0;
		this->IN.load[0].F[1] = 0;
		this->IN.load[0].F[2] = 0;
		this->IN.load[0].T[0] = 0;
		this->IN.load[0].T[1] = 0;
		this->IN.load[0].T[2] = 0;

		// 潤滑理論
		this->IN.rigid.g[0] = 0;
		this->IN.rigid.g[1] = 0;
		this->IN.rigid.g[2] = 0;
		this->IN.rigid.l = 1;
		this->IN.rigid.t = 1;
		this->IN.bound.v_const[0] = false;
		this->IN.bound.v_const[1] = false;
		this->IN.bound.v_const[2] = false;
		this->IN.bound.w_const[0] = false;
		this->IN.bound.w_const[1] = false;
		this->IN.bound.w_const[2] = false;
		this->BS.allocate(IN);
		this->OUT.allocate(this->IN);
		double v0 = 0, w0 = 1, wn = 0;
		this->BS.init(this->IN, v0, w0, wn);

	}

	virtual void TearDown() {
	}
};

// (玉初期値計算) 玉が等間隔に並んでいるか確認
TEST_F(BS_BallSrewTest, init_position_case) {
	this->SingleNut_init();
	this->BS.LD.F = Vector3d(100, 0, 0);
	this->BS.save(this->OUT);
	double l = 181.2879e-3;
	for (int i = 0; i < 5; i++) {
		double x = this->OUT.CC[0].BL[i].x[0];
		double y = this->OUT.CC[0].BL[i].x[1];
		double z = this->OUT.CC[0].BL[i].x[2];
		double _x = l / 4 * i;
		EXPECT_NEAR(x, _x, 1e-6);
		double _th = Numeric::pi * 2 / 4 * i;
		double _y = cos(_th) * 50e-3;
		double _z = sin(_th) * 50e-3;
		EXPECT_NEAR(y, _y, 1e-6);
		EXPECT_NEAR(z, _z, 1e-6);
	}
	return;
}

// (初期値計算ケース0) 荷重0のときは何もしない．
TEST_F(BS_BallSrewTest, preset_y0_case0) {
	this->SingleNut_init();
	this->BS.LD.F = Vector3d(0, 0, 0);
	this->BS.preset_y0(1e-9, 1e-9);
	this->BS.save(this->OUT);

	double x = this->OUT.ST.x[0];
	double y = this->OUT.ST.x[1];
	double z = this->OUT.ST.x[2];
	EXPECT_NEAR(x, 0, 1e-6);
	EXPECT_NEAR(y, 0, 1e-6);
	EXPECT_NEAR(z, 0, 1e-6);
	return;
}
// (初期値計算ケース1) 純アキシアルのとき，大体アキシアルすきま分シャフトと玉が移動するか確認
TEST_F(BS_BallSrewTest, preset_y0_case1) {
	this->SingleNut_init();
	this->BS.LD.F = Vector3d(100, 0, 0);
	this->BS.preset_y0(1e-9, 1e-9);
	this->BS.save(this->OUT);

	double x = this->OUT.ST.x[0];
	double y = this->OUT.ST.x[1];
	double z = this->OUT.ST.x[2];
	EXPECT_NEAR(x, 0.86602540e-3 * 2, 0.2e-3);
	EXPECT_NEAR(y, 0, 1e-6);
	EXPECT_NEAR(z, 0, 1e-6);

	double l = 181.2879e-3;		// リード長[m]
	for (int i = 0; i < 5; i++) {
		double x = this->OUT.CC[0].BL[i].x[0];
		double y = this->OUT.CC[0].BL[i].x[1];
		double z = this->OUT.CC[0].BL[i].x[2];
		double _x = l / 4 * i + 0.86602540e-3;
		double _th = Numeric::pi * 2 / 4 * i;
		double _y = cos(_th) * 50.5e-3;
		double _z = sin(_th) * 50.5e-3;
		// 許容誤差はy, z座標の1%程度
		EXPECT_NEAR(this->OUT.CC[0].BL[i].x[0], _x, 5e-4);
		EXPECT_NEAR(this->OUT.CC[0].BL[i].x[1], _y, 5e-4);
		EXPECT_NEAR(this->OUT.CC[0].BL[i].x[2], _z, 5e-4);
	}
	return;
}

// (初期値計算ケース2) 純ラジアルのとき，大体ラジアルすきま分シャフトが移動するか確認
// また，すべての玉がナットに接触しており，等間隔になっているか確認
TEST_F(BS_BallSrewTest, preset_y0_case2) {
	this->SingleNut_init();
	this->BS.LD.F = Vector3d(0, 100, 0);
	this->BS.preset_y0(1e-9, 1e-9);
	this->BS.save(this->OUT);

	double x = this->OUT.ST.x[0];
	double y = this->OUT.ST.x[1];
	double z = this->OUT.ST.x[2];
	EXPECT_NEAR(x, 0, 1e-6);
	EXPECT_NEAR(y, 1e-3, 0.1e-3);
	EXPECT_NEAR(z, 0, 1e-6);

	double l = 181.2879e-3;		// リード長[m]
	for (int i = 0; i < 5; i++) {
		double x = this->OUT.CC[0].BL[i].x[0];
		double y = this->OUT.CC[0].BL[i].x[1];
		double z = this->OUT.CC[0].BL[i].x[2];
		double _x = l / 4 * i;
		double _th = Numeric::pi * 2 / 4 * i;
		double _y = cos(_th) * 50.5e-3;
		double _z = sin(_th) * 50.5e-3;
		EXPECT_NEAR(x, _x, 5e-4);
		EXPECT_NEAR(y, _y, 5e-4);
		EXPECT_NEAR(z, _z, 5e-4);
	}
	return;
}

// (初期値計算ケース3) ミスアライメントがあるとき，荷重をかけた向きに軸が傾くか確認
TEST_F(BS_BallSrewTest, preset_y0_case3) {
	this->SingleNut_init();
	this->BS.LD.F = Vector3d(0, 0, 0);
	this->BS.LD.T = Vector3d(0, 10, 0);
	this->BS.preset_y0(1e-9, 1e-9);
	this->BS.save(this->OUT);

	double x = this->OUT.ST.x[0];
	double y = this->OUT.ST.x[1];
	double z = this->OUT.ST.x[2];
	double ax = this->OUT.ST.ax[0];
	double ay = this->OUT.ST.ax[1];
	double az = this->OUT.ST.ax[2];
	EXPECT_NEAR(ay, 0, 1e-6);
	EXPECT_LT(az, 0);
	return;
}