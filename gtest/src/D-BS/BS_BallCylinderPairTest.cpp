/*******************************************************************************
!						"BS_BallCylinderPairTest.cpp"
!														2021/05/13	[Core-T] 柏原
!	基本的に D-B4P の B4P_BallRingPairTest.cpp のテストケースを流用
!
!*******************************************************************************/

#include "pch.h"
#include "BS_FileIn_stab.h"
#include "BS_FileOut_stab.h"
using Eigen::IOFormat;
using namespace std;

#define cos30 0.86602540378443864676372317075294
#define sin30 0.5

class BS_BallCylinderPairTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	BS_BallCylinderPair BSP, BNP;
	BS_FileIn_stab IN;
	BS_FileOut_stab FO;
	BS_Shaft ST;
	BS_SingleNut NT;
	Ball *BL;

	virtual void SetUp() {
	}
	void init() {
		this->IN.shaft.resize(1);
		this->IN.shaft[0].density = 1.0;
		this->IN.shaft[0].poisson = 1.0;
		this->IN.shaft[0].young = 1.0;
		this->IN.shaft[0].ri = 1.2;
		this->IN.shaft[0].ro = 1.2;

		this->IN.shaft[0].spiral.resize(1);

		this->IN.shaft[0].spiral[0].alp = 0.0;
		this->IN.shaft[0].spiral[0].l = 0.8;
		this->IN.shaft[0].spiral[0].r = 1.0;

		for (int i = 0; i < 2; i++) {
			this->IN.shaft[0].spiral[0].groove[i].eta[0] = 0.0;
			this->IN.shaft[0].spiral[0].groove[i].eta[1] = 0.0;
			this->IN.shaft[0].spiral[0].groove[i].r = 1.0e-3;
			this->IN.shaft[0].spiral[0].groove[i].sigma = 1.0e-6;
		}

		this->ST.allocate(this->IN.shaft);

		double v0 = 0.0;
		double w0 = 0.0;
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}
		this->ST.init(this->IN.shaft, x_const, Rx_const, v0, w0);

		Rigid::l = 1e-2;
		Rigid::t = 1e-3;

		this->ST.init_pos(v0, w0);

		this->ST.x = Vector3d(0.0, 1.0, 0.0);
		this->ST.set_dx();
		this->BSP.init(IN.BallShaftPair[0], IN.tribology, IN.oil);
	}

	// ボールねじの幾何形状を定義
	// 溝一つ，螺旋の傾きが30°の単純形状にする
	void BS_init() {
		// Pairクラス物性値
		this->IN.BallShaftPair.resize(1);
		this->IN.BallShaftPair[0].groove[0].zeta = 0.2;	// 減衰項
		this->IN.BallShaftPair[0].groove[0].mu = 0.1;		// 摩擦係数
		this->IN.BallShaftPair[0].groove[1].zeta = 0.2;	// 減衰項
		this->IN.BallShaftPair[0].groove[1].mu = 0.1;		// 摩擦係数
		// 玉物性値
		double D = 20e-3;				// 玉径[m]
		double E = 207760000000;		// ヤング率[Pa]
		double por = 0.3;				// ポアソン比[-]
		double den = 7830;				// 密度[kg/m^3]
		double rms = 0.00000002;		// 粗さrms[m]
		this->BL = new Ball;
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}
		this->BL->init(D, E, por, den, rms, x_const, Rx_const);

		// シャフト物性値
		IN.shaft.resize(1);
		IN.shaft[0].spiral.resize(1);
		IN.shaft[0].spiral[0].groove[0].eta[1] = 0;
		IN.shaft[0].spiral[0].groove[1].eta[1] = 0;
		IN.shaft[0].spiral[0].groove[0].r = 11e-3;
		IN.shaft[0].spiral[0].groove[1].r = 1;
		// 溝R中心PCD[m]（溝pcd半径，0.0357535 m くらい）
		IN.shaft[0].spiral[0].groove[0].eta[0] = 0;
		IN.shaft[0].spiral[0].groove[1].eta[0] = 0;

		IN.shaft[0].spiral[0].groove[0].sigma = 0.00000006;
		IN.shaft[0].spiral[0].groove[1].sigma = 0.00000006;
		IN.shaft[0].spiral[0].alp = 0.5 * Numeric::pi;
		IN.shaft[0].spiral[0].r = 50e-3;
		IN.shaft[0].spiral[0].l = 181.2879e-3;// (= 2πr * tan30°)

		IN.shaft[0].m = 1;
		IN.shaft[0].young = 207760000000;
		IN.shaft[0].poisson = 0.3;
		IN.shaft[0].Ix = 1, IN.shaft[0].Iyz = 1, IN.shaft[0].Iyz = 1;

		this->IN.oil.eta = 0.5;	// 粘度（Pa*s）(適当)
		this->IN.oil.beta = 0.5;	// 温度粘度係数（適当）
		this->IN.oil.k = 0.145;	// 油熱伝導率[W/(m・K)]
		this->IN.oil.alpha = 0.2;	// 圧力粘度係数[mm2/kgf]
		this->IN.oil.lm = -1;	// メニスカス長さ


		// 潤滑理論
		IN.tribology.rollingresistance = BS_In::Tribology::RollingResistance::RollingResistanceNothing;
		IN.tribology.coulomb = BS_In::Tribology::Tangent;
		IN.tribology.coulomb_slope = 1000;
		IN.tribology.filmThickness = BS_In::Tribology::FilmThickness::FilmThicknessNothing;
		IN.tribology.hysteresis = BS_In::Tribology::HysteresisNothing;
		IN.tribology.hysteresis_factor = 0;

		double v0 = 0.0;
		double w0 = 0.0;
		this->ST.CY = new BS_Cylinder[1];
		this->ST.CY[0].allocate(1);

		this->ST.init(this->IN.shaft, x_const, Rx_const, v0, w0);

		Rigid::l = 1;
		Rigid::t = 1;

		this->ST.init_pos(v0, w0);
		this->BSP.link(this->BL, this->ST.CY, 0);
		this->BSP.init(IN.BallShaftPair[0], IN.tribology, IN.oil);


		// 外輪物性値
		this->IN.nut.resize(1);
		this->IN.nut[0].spiral.resize(1);

		this->IN.nut[0].spiral[0].groove[0].eta[1] = 0;
		this->IN.nut[0].spiral[0].groove[1].eta[1] = 0;
		this->IN.nut[0].spiral[0].groove[0].r = 11e-3;
		this->IN.nut[0].spiral[0].groove[1].r = 1e3;
		this->IN.nut[0].spiral[0].groove[0].eta[0] = 0;
		this->IN.nut[0].spiral[0].groove[1].eta[0] = 0;
		this->IN.nut[0].spiral[0].groove[0].sigma = 0.00000006;
		this->IN.nut[0].spiral[0].groove[1].sigma = 0.00000006;
		IN.nut[0].spiral[0].alp = 0.5 * Numeric::pi;
		this->IN.nut[0].spiral[0].r = 50e-3;
		IN.nut[0].spiral[0].l = 181.2879e-3;// (= 2πr * tan30°)

		this->IN.nut[0].m = 1;
		this->IN.nut[0].young = 207760000000;
		this->IN.nut[0].poisson = 0.3;
		this->IN.nut[0].Ix = 1, this->IN.nut[0].Iyz = 1, this->IN.nut[0].Iyz = 1;

		this->NT.CY = new BS_Cylinder[1];
		this->NT.CY[0].allocate(1);
		this->NT.init(this->IN.nut, w0);
		this->BNP.link(this->BL, this->NT.CY, 0);
		this->BNP.init(IN.BallShaftPair[0], IN.tribology, IN.oil);
		this->IN.circuit.resize(1);
		this->IN.circuit[0].ball.resize(1);

		this->FO.allocate(this->IN);
		this->FO.BNP.resize(1);
		//this->BNP.TB.TR = new Tribology::AiharaT();
		//this->BNP.TB.CL = new Tribology::Tangent();
		//this->BNP.TB.RR = new Tribology::RollingResistanceNothing();
		//this->BNP.TB.FT = new Tribology::FilmThicknessNothing();	// 油膜厚さを0に固定
		//this->BNP.TB.cs = 1000;
	}

	// 4玉の軸受緒元（25BSWZ01）を流用
	void init_25BSWZ() {
		// 入力データは下記のものを用いた
		// \\ans00978\kiken\BRAIN\ソースコード他\DBRAINpp\docs\06_D4Bv100\05_単体テスト\B4P_BallRingPair
		// 軸受型番　25BSWZ01
		// Pairクラス物性値
		this->IN.BallShaftPair.resize(1);
		this->IN.BallShaftPair[0].groove[0].zeta = 0.2;	// 減衰項
		this->IN.BallShaftPair[0].groove[0].mu = 0.1;		// 摩擦係数
		this->IN.BallShaftPair[0].groove[1].zeta = 0.2;	// 減衰項
		this->IN.BallShaftPair[0].groove[1].mu = 0.1;		// 摩擦係数
		// 玉物性値
		double D = 0.00635;				// 玉径[m]
		double E = 207760000000;		// ヤング率[Pa]
		double por = 0.29;				// ポアソン比[-]
		double den = 7830;				// 密度[kg/m^3]
		double rms = 0.00000002;		// 粗さrms[m]
		this->BL = new Ball;
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}
		this->BL->init(D, E, por, den, rms, x_const, Rx_const);

		// 内輪物性値
		IN.shaft.resize(1);
		IN.shaft[0].spiral.resize(1);
		IN.shaft[0].spiral[0].r = 0.0355 / 2;
		IN.shaft[0].spiral[0].groove[0].eta[1] = 0.0000725;
		IN.shaft[0].spiral[0].groove[1].eta[1] = -0.0000725;
		IN.shaft[0].spiral[0].groove[0].r = 0.00332105;
		IN.shaft[0].spiral[0].groove[1].r = 0.00332105;
		// 溝R中心PCD[m]（溝pcd半径，0.0357535 m くらい）
		IN.shaft[0].spiral[0].groove[0].eta[0] = -(0.035753569339629220 - 0.0355) / 2;
		IN.shaft[0].spiral[0].groove[1].eta[0] = -(0.035753569339629220 - 0.0355) / 2;

		IN.shaft[0].spiral[0].groove[0].sigma = 0.00000006;
		IN.shaft[0].spiral[0].groove[1].sigma = 0.00000006;
		IN.shaft[0].spiral[0].alp = 0.5 * Numeric::pi;
		IN.shaft[0].spiral[0].l = 0;

		IN.shaft[0].m = 1;
		IN.shaft[0].young = 207760000000;
		IN.shaft[0].poisson = 0.29;
		IN.shaft[0].Ix = 1, IN.shaft[0].Iyz = 1, IN.shaft[0].Iyz = 1;

		this->IN.oil.eta = 0.5;	// 粘度（Pa*s）(適当)
		this->IN.oil.beta = 0.5;	// 温度粘度係数（適当）
		this->IN.oil.k = 0.145;	// 油熱伝導率[W/(m・K)]
		this->IN.oil.alpha = 0.2;	// 圧力粘度係数[mm2/kgf]
		this->IN.oil.lm = -1;	// メニスカス長さ


		// 潤滑理論
		IN.tribology.rollingresistance = BS_In::Tribology::RollingResistance::RollingResistanceNothing;
		IN.tribology.coulomb = BS_In::Tribology::Tangent;
		IN.tribology.coulomb_slope = 1000;
		IN.tribology.filmThickness = BS_In::Tribology::FilmThickness::FilmThicknessNothing;
		IN.tribology.hysteresis = BS_In::Tribology::HysteresisNothing;
		IN.tribology.hysteresis_factor = 0;

		double v0 = 0.0;
		double w0 = 0.0;
		this->ST.CY = new BS_Cylinder[1];
		this->ST.CY[0].allocate(1);
		this->ST.init(this->IN.shaft, x_const, Rx_const, v0, w0);

		Rigid::l = 1;
		Rigid::t = 1;

		this->ST.init_pos(v0, w0);
		this->BSP.link(this->BL, this->ST.CY, 0);
		this->BSP.init(IN.BallShaftPair[0], IN.tribology, IN.oil);


		// 外輪物性値
		this->IN.nut.resize(1);
		this->IN.nut[0].spiral.resize(1);
		this->IN.nut[0].spiral[0].r = 0.0355 / 2;

		this->IN.nut[0].spiral[0].groove[0].eta[1] = 0.0000725;
		this->IN.nut[0].spiral[0].groove[1].eta[1] = -0.0000725;
		this->IN.nut[0].spiral[0].groove[0].r = 0.00332105;
		this->IN.nut[0].spiral[0].groove[1].r = 0.00332105;
		this->IN.nut[0].spiral[0].groove[0].eta[0] = (0.035246430660370774 - 0.0355) / 2;;
		this->IN.nut[0].spiral[0].groove[1].eta[0] = (0.035246430660370774 - 0.0355) / 2;;
		this->IN.nut[0].spiral[0].groove[0].sigma = 0.00000006;
		this->IN.nut[0].spiral[0].groove[1].sigma = 0.00000006;
		IN.nut[0].spiral[0].alp = 0.5 * Numeric::pi;
		IN.nut[0].spiral[0].l = 0.1;

		this->IN.nut[0].m = 1;
		this->IN.nut[0].young = 207760000000;
		this->IN.nut[0].poisson = 0.29;
		this->IN.nut[0].Ix = 1, this->IN.nut[0].Iyz = 1, this->IN.nut[0].Iyz = 1;

		this->NT.CY = new BS_Cylinder[1];
		this->NT.CY[0].allocate(1);
		this->NT.init(this->IN.nut, w0);
		this->BNP.link(this->BL, this->NT.CY, 0);
		this->BNP.init(IN.BallShaftPair[0], IN.tribology, IN.oil);
		this->FO.BNP.resize(1);

		this->IN.circuit.resize(1);
		this->IN.circuit[0].ball.resize(1);
		this->FO.allocate(this->IN);
	}


	virtual void TearDown() {
	}
};

// (荷重計算ケース1) シンプルなケースを想定
// 接触角180°かつ玉位相角0°で無回転で壁面に対して並進運動している時を仮定
TEST_F(BS_BallCylinderPairTest, calc_force_case1) {
	this->BS_init();
	// 玉接近量 1um，玉の進行方向は螺旋接線方向になるように設定
	Vector3d bl_x = Vector3d(0, 0, 50e-3 + 1.001e-3);
	Vector3d bl_v = Vector3d(0.5, 0.8660254, 0);
	Vector3d bl_w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）

	this->BL->set_y(bl_x, bl_v, q, bl_w);


	Vector3d Fbi, Tbi, Fib, Tib, Fs, Ts;

	this->BNP.get_FT(Fbi, Tbi, Tib, Fs, Ts);

	// 1. 玉からナットに働くヘルツ荷重は -η 方向（= +z 方向）の成分だけになっているか確認
	Vector3d Fn = -Fbi - Fs;
	double err = 1e-3;
	EXPECT_NEAR(Fn.z(), Fn.norm(), Fn.norm() * err);
	// 2. 摩擦力 = ヘルツ荷重 × 摩擦係数 になっているか確認
	double F = Fs.norm();
	EXPECT_NEAR(Fs.norm(), Fn.norm() * 0.1, Fs.norm() * err);
	// 3. 玉からナットに働く摩擦は玉が進行する向きとは同じ向きになっていることを確認
	EXPECT_NEAR((Fs.normalized() - bl_v.normalized()).norm(), 0, err);


	return;
}

// (荷重計算ケース2) 位相角・接触角が非ゼロのケースを検証
// 接触角150°かつ玉位相角90°で無回転で壁面に対して並進運動している時を仮定
TEST_F(BS_BallCylinderPairTest, calc_force_case2) {
	this->BS_init();
	// 玉接近量 1um，玉の進行方向は螺旋接線方向になるように設定
	// この螺旋の初期位相角は90°なので，軸から見て-y方向に玉が存在
	Vector3d bl_x = Vector3d(181.2879e-3 * 0.25 + 1.001e-3 * sin30 * cos30, -50e-3 - 1.001e-3 * cos30, 1.001e-3 * sin30 * sin30);
	Vector3d bl_v = Vector3d(sin30, 0, -cos30);
	Vector3d bl_w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）

	this->BL->set_y(bl_x, bl_v, q, bl_w);


	Vector3d Fbi, Tbi, Fib, Tib, Fs, Ts;

	this->BNP.get_FT(Fbi, Tbi, Tib, Fs, Ts);

	// 1. 玉からナットに働くヘルツ荷重はηに対して150°の向きになっているか確認
	Vector3d e = (-Fbi - Fs).normalized();
	double err = 1e-3;
	Vector3d e_ex = Vector3d(sin30 * cos30, -cos30, sin30 * sin30);
	IOFormat CleanFmt(4, 0, ", ", ",", "", "", "[", "]");
	EXPECT_NEAR((e_ex - e).norm(), 0, err) << "e = " << e.format(CleanFmt) << endl << "e_ex = " << e_ex.format(CleanFmt);
	// 2. 摩擦力 = ヘルツ荷重 × 摩擦係数 になっているか確認
	Vector3d Fn = -Fbi - Fs;
	EXPECT_NEAR(Fs.norm(), Fn.norm() * 0.1, Fs.norm() * err);
	// 3. 玉からナットに働く摩擦は玉が進行する向きとは同じ向きになっていることを確認
	EXPECT_NEAR((Fs.normalized() - bl_v.normalized()).norm(), 0, err);

	return;
}

// (荷重計算ケース3) 玉がスピン方向に回転しているケースを検証
// 接触角0°かつ玉位相角0°で玉が壁面に対してスピンしている時を仮定
TEST_F(BS_BallCylinderPairTest, calc_force_case3) {
	this->BS_init();
	// 玉接近量 1um，玉の進行方向は螺旋接線方向になるように設定
	// この螺旋の初期位相角は90°なので，軸から見て+z方向に玉が存在
	Vector3d bl_x = Vector3d(0, 0, 50e-3 + 1.001e-3);
	Vector3d bl_v = Vector3d(0, 0, 0);
	Vector3d bl_w = Vector3d(0, 0, 1);		// -η 方向に自転
	Quaterniond q = Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）

	this->BL->set_y(bl_x, bl_v, q, bl_w);
	Vector3d Fbi, Tbi, Fib, Tib, Fs, Ts;

	// 荷重計算し，スライスの結果を取得
	this->BNP.get_FT(Fbi, Tbi, Tib, Fs, Ts);
	this->BNP.save(FO.BNP[0]);

	// 1. 各スライス片のすべり摩擦の正負を確認
	// 接触楕円の -ζ 側にあるスライス片のすべり摩擦は -ξ 方向
	for (int i = 0; i < 11; i++) {
		double fs_xai = FO.BNP[0].GV[0].SL[i].fs_[0];
		double ps_zeta = FO.BNP[0].GV[0].SL[i].ps_[2];
		EXPECT_LT(fs_xai, 0);
	}
	// 接触楕円の +ζ 側にあるスライス片のすべり摩擦は +ξ 方向
	for (int i = 12; i < 21; i++) {
		double fs_xai = FO.BNP[0].GV[0].SL[i].fs_[0];
		double ps_zeta = FO.BNP[0].GV[0].SL[i].ps_[2];
		EXPECT_GT(fs_xai, 0);
	}
	// 2. 摩擦力 << ヘルツ荷重 × 摩擦係数 になっているか確認
	Vector3d Fn = -Fbi - Fs;
	EXPECT_LT(Fs.norm(), Fn.norm() * 0.1);

	// 3. スライス10 の影響を差し引けば，玉自転とは逆向きに摩擦が発生することを確認
	Vector3d Fs__, Ts__, ts_, fs_;
	// スライス10の摩擦は螺旋座標系で出力されるので，慣性座標系に直してから出力
	Matrix3d xyz2eta = this->BNP.CY->get_xyz2eta(BNP.iSP, this->BNP.SV.eta[0]);
	for (int i = 0; i < 3; i++) {
		ts_[i] = FO.BNP[0].GV[0].SL[10].ts_[i];
		fs_[i] = FO.BNP[0].GV[0].SL[10].fs_[i];
	}
	Vector3d ts = this->BNP.CY->to_inertialvector(ts_, xyz2eta);
	Vector3d fs = this->BNP.CY->to_inertialvector(fs_, xyz2eta);
	Ts__ = Tbi - ts;
	Fs__ = -Fs - fs;

	EXPECT_LT(Fs__.norm(), 1e-2);
	EXPECT_LT((Ts__.normalized() + bl_w.normalized()).norm(), 1e-2);

	return;
}


// (静解析step0) calc_force_case1 と同じ条件で荷重計算し，結果を確認
TEST_F(BS_BallCylinderPairTest, get_F0_case1) {
	this->BS_init();
	// 玉接近量 1um，玉の進行方向は螺旋接線方向になるように設定
	Vector3d bl_x = Vector3d(0, 0, 50e-3 + 1.001e-3);
	Vector3d bl_v = Vector3d(0.5, 0.8660254, 0);
	Vector3d bl_w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）

	this->BL->set_y(bl_x, bl_v, q, bl_w);

	Vector2d Fbc;
	Vector3d Fcb, Tcb;
	this->BNP.get_F0(Fbc, Fcb, Tcb);

	// 1. 玉からナットに働くヘルツ荷重は -η 方向（= +z 方向）の成分だけになっているか確認
	double err = 1e-3;
	EXPECT_NEAR(Fbc[0], Fbc.norm(), Fbc.norm() * err);
	EXPECT_NEAR(Fcb[2], Fcb.norm(), Fcb.norm() * err);
	return;
}

// (静解析step2) calc_force_case1 と同じ条件で荷重計算し，結果を確認
TEST_F(BS_BallCylinderPairTest, get_F1_case1) {
	this->BS_init();
	// 玉接近量 1um，玉の進行方向は螺旋接線方向になるように設定
	Vector3d bl_x = Vector3d(0, 0, 50e-3 + 1.001e-3);
	Vector3d bl_v = Vector3d(0.5, 0.8660254, 0);
	Vector3d bl_w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）

	this->BL->set_y(bl_x, bl_v, q, bl_w);

	Vector2d Fbc;
	Vector3d Fcb, Tcb;
	this->BNP.get_F1(true, Fbc, Fcb, Tcb);

	// ステップ1 の接触荷重はステップ0の荷重に摩擦荷重を加えたものなので，
	// ステップ0の結果に予想される摩擦荷重を加算してステップ1の予測値を導出．
	Vector3d Fcb_ex = Vector3d(-3.1495739376147817 * cos30, -3.1495739376147817 * sin30, 31.495739376147817);
	EXPECT_LT((Fcb - Fcb_ex).norm(), 1e-2);
	return;
}
// (静解析step2) シャフト側で荷重計算し，結果を確認
TEST_F(BS_BallCylinderPairTest, get_F1_case2) {
	this->BS_init();
	// 玉接近量 1um，玉の進行方向は螺旋接線方向になるように設定
	Vector3d bl_x = Vector3d(0, 0, 50e-3 - 1.001e-3);
	Vector3d bl_v = Vector3d(0.5, 0.8660254, 0);
	Vector3d bl_w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）

	this->BL->set_y(bl_x, bl_v, q, bl_w);
	Vector2d Fbc;
	Vector3d Fcb, Tcb;
	this->BSP.get_F1(false, Fbc, Fcb, Tcb);

	// ステップ1 の接触荷重はステップ0の荷重に摩擦荷重を加えたものなので，
	// ステップ0の結果に予想される摩擦荷重を加算してステップ1の予測値を導出．
	Vector3d Fcb_ex = Vector3d(-3.0082504011194946 * cos30, -3.0082504011194946 * sin30, -30.082504011194946);
	EXPECT_LT((Fcb - Fcb_ex).norm(), 1e-2);
	return;
}

// (静解析step3) calc_force_case1 と同じ条件で荷重計算し，結果を確認
TEST_F(BS_BallCylinderPairTest, get_F2_case1) {
	this->BS_init();
	// 玉接近量 1um，玉の進行方向は螺旋接線方向になるように設定
	Vector3d bl_x = Vector3d(0, 0, 50e-3 + 1.001e-3);
	Vector3d bl_v = Vector3d(sin30, cos30, 0);
	Vector3d bl_w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）

	this->BL->set_y(bl_x, bl_v, q, bl_w);
	Vector2d Fbc;
	Vector3d Fcb, Tcb;
	this->BNP.get_F1(true, Fbc, Fcb, Tcb);

	Vector3d vF, vT;
	this->BNP.get_F2(vF, vT);

	// ステップ2 の接触荷重はステップ0の荷重にすべり速度をかけたものになっているか確認
	Vector3d vF_ex = Vector3d(31.495739376147817 * sin30, 31.495739376147817 * cos30, 0);
	EXPECT_LT((vF - vF_ex).norm(), 1e-2);
	return;
}