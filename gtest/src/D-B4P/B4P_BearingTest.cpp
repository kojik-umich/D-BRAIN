#include "pch.h"
#include "B4P_StabPair.h"
#include "B4P_StabIn.h"
#include "B4P_StabOut.h"
using Eigen::Matrix2d;


class B4P_BearingTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	B4P_Bearing b4p;
	B4P_StabIn FI;
	B4P_StabOut FO;

protected:
	virtual void SetUp() {


		return;
	}

	virtual void TearDown() {

	}
	// テスト用の初期化関数
	void test_init(int Z) {
		this->b4p.Z = Z;
		int msmax = 21;
		this->FO.allocate(Z, _MAX_CONTACT_, msmax);
		this->b4p.BOP = new B4P_StabBallOuterRingPair[Z];
		this->b4p.BIP = new B4P_StabBallInnerRingPair[Z];
		this->b4p.BCP = new B4P_StabBallCagePair[Z];
		this->b4p.BL = new Ball[Z];
		this->b4p.CG = new bal_SnapCage();

		double D = 0.00635;				// 玉径[m]
		double E = 207760000000;		// ヤング率[Pa]
		double por = 0.29;				// ポアソン比[-]
		double den = 7830;				// 密度[kg/m^3]
		double rms = 0.00000002;		// 粗さrms[m]
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}
		this->b4p.BL[0].init(D, E, por, den, rms, x_const, Rx_const);
		this->b4p.BL[0].m = 1.0;
		this->b4p.IR.m = 1.0;
		this->b4p.IR.m_inv = 1.0;
		this->b4p.IR.I = Vector3d::Ones();
		this->b4p.IR.I_inv = Vector3d::Ones();
		this->b4p.CG->m = 1.0;
		this->b4p.CG->m_inv = 1.0;
		this->b4p.CG->I = Vector3d::Ones();
		this->b4p.CG->I_inv = Vector3d::Ones();

		this->b4p.F_load = Vector3d::Ones();
		this->b4p.T_load = Vector3d::Ones();
		FI.BL = new  B4P_In::Ball[1];
		// 入力データは下記のものを用いた
		// \\ans00978\kiken\BRAIN\ソースコード他\DBRAINpp\docs\06_D4Bv100\05_単体テスト\B4P_BallRingPair
		// 軸受型番　25BSWZ01
		// Pairクラス物性値
		FI.BOP.dzeta = 0.2;	// 減衰項
		FI.BOP.mu = 0.1;		// 摩擦係数
		FI.BIP.dzeta = 0.2;	// 減衰項
		FI.BIP.mu = 0.1;		// 摩擦係数

		// 玉物性値
		FI.balldia = 0.00635;				// 玉径[m]
		FI.ballnum = 1;
		FI.BL[0].dia = 0.00635, FI.BL[0].E = 207760000000, FI.BL[0].por = 0.29, FI.BL[0].den = 7830, FI.BL[0].rms = 0.00000002;
		// 内輪物性値
		FI.ballpcd = 0.0355;
		double clr = 0;

		FI.IR.Rox[0] = -0.0000725;
		FI.IR.Rox[1] = 0.0000725;
		FI.IR.R[0] = 0.00332105;
		FI.IR.R[1] = 0.00332105;
		FI.IR.Rod[0] = FI.ballpcd - clr * 0.5 + 2 * sqrt((FI.IR.R[0] - D * 0.5)*(FI.IR.R[0] - D * 0.5) - FI.IR.Rox[0] * FI.IR.Rox[0]);
		FI.IR.Rod[1] = FI.ballpcd - clr * 0.5 + 2 * sqrt((FI.IR.R[1] - D * 0.5)*(FI.IR.R[1] - D * 0.5) - FI.IR.Rox[1] * FI.IR.Rox[1]);
		FI.IR.hedge[0];
		FI.IR.hedge[1];
		FI.IR.rms = 0.00000006;
		FI.IR.m = 1;
		FI.IR.E = 207760000000;
		FI.IR.por = 0.29;
		FI.IR.Ix, FI.IR.Iyz, FI.IR.Iyz;

		// 潤滑理論
		FI.TB.rollingresistance;

		// 外輪物性値
		FI.OR.Rox[0] = -0.0000725;
		FI.OR.Rox[1] = 0.0000725;
		FI.OR.R[0] = 0.00332105;
		FI.OR.R[1] = 0.00332105;
		FI.OR.Rod[0] = FI.ballpcd + clr * 0.5 - 2 * sqrt((FI.OR.R[0] - D * 0.5)*(FI.OR.R[0] - D * 0.5) - FI.OR.Rox[0] * FI.OR.Rox[0]);
		FI.OR.Rod[1] = FI.ballpcd + clr * 0.5 - 2 * sqrt((FI.OR.R[1] - D * 0.5)*(FI.OR.R[1] - D * 0.5) - FI.OR.Rox[1] * FI.OR.Rox[1]);
		FI.OR.hedge[0];
		FI.OR.hedge[1];
		FI.OR.rms = 0.00000006;
		FI.OR.m = 1.0;
		FI.OR.E = 2.0776e11;
		FI.OR.por = 0.29;
		FI.OR.Ix, FI.OR.Iyz, FI.OR.Iyz;
		FI.LB.eta0 = 0.5;	// 粘度（Pa*s）(適当)
		FI.LB.beta0 = 0.5;	// 温度粘度係数（適当）
		FI.LB.k0 = 0.145;	// 油熱伝導率[W/(m・K)]
		FI.LB.alpha0 = 0.2;	// 圧力粘度係数[mm2/kgf]
		FI.LB.lm0 = -1;		// メニスカス長さ

		Vector3d x = Vector3d(0, 0, 0);
		Vector3d v = Vector3d(0, 0, 0);
		Quaterniond q = Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）
		Vector3d w = Vector3d(0, 0, 0);

		Rigid::l = 1.0;
		Rigid::t = 1.0;
		Rigid::g = Vector3d(0, 0, 9.8);
		FI.TB.filmThickness = B4P_In::Tribology::HamrockDowsonHc;

		FI.TB.coulomb = B4P_In::Tribology::Tangent;
		FI.TB.coulomb_slope = 1e3;
		FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::Aihara;
		FI.TB.hysteresis = B4P_In::Tribology::HysteresisNothing;
		double x0[3], ax0[3];
		for (int i = 0; i < 3; i++)
			x0[i] = 0.0;
		ax0[0] = 1.0; ax0[1] = 0.0; ax0[2] = 0.0;

		this->b4p.init(FI,  x0, ax0);
		return;
	}
};

TEST_F(B4P_BearingTest, get_dydt0) {
	this->test_init(1);

	this->b4p.BL[0].x = Vector3d(0, 0, 0);
	this->b4p.CG->x = Vector3d(0, 0, 0);
	this->b4p.IR.x = Vector3d(0, 0, 0);
	this->b4p.OR.x = Vector3d(0, 0, 0);
	double dydt[52];
	this->b4p.get_dydt(dydt);
	Vector3d bl_dydt_ex = Vector3d(0, 0, 9.8);
	Vector3d bl_dydt = Vector3d(dydt[26], dydt[27], dydt[28]);

	ASSERT_TRUE(false);
	return;
};

// 前半：initのテスト，後半：make_AllParamsのテスト
TEST_F(B4P_BearingTest, init) {

	double dydt[52];

	B4P_StabIn FI;

	FI.cos_alp0 = cos(30 * 0.017453292519943295);
	FI.balldia = 0.00635, FI.ballnum = 12, FI.ballpcd = 0.0355, FI.omegair = 31.4, FI.omegaor = 0;
	bool v_is_Locked = true;
	bool w_is_Locked = true;
	FI.cage_type = B4P_In::snap_cage;
	FI.Cage.rmg0[0] = 0;
	FI.Cage.rmg0[1] = 0;
	FI.Cage.rmg0[2] = 0;
	FI.LoadIn[0] = 100, FI.LoadIn[1] = 0, FI.LoadIn[2] = 200;
	FI.LoadIn[3] = 0, FI.LoadIn[4] = 60, FI.LoadIn[5] = 0;
	FI.rigid.g[0] = 0, FI.rigid.g[1] = 0, FI.rigid.g[2] = 9.8;
	double x0[3], ax0[3];
	x0[0] = 1.0; x0[1] = 0.0; x0[2] = 0.0;
	ax0[0] = 1.0; ax0[1] = 0.0; ax0[2] = 1.0;

	this->test_init(FI.ballnum);

	FI.BL = new B4P_In::Ball[FI.ballnum];
	for (int i = 0; i < FI.ballnum; i++) {
		FI.BL[i].den = 7830;
		FI.BL[i].E = 208000000000;
		FI.BL[i].por = 0.29;
		FI.BL[i].rms = 0.00000006;
		FI.BL[i].dia = 0.00635;
	}
	this->b4p.init(FI, x0, ax0);
	double err = 0.01;		// 許容誤差

	// 玉位置の検証
	double half_pcd = 0.0355 * 0.5;
	Vector3d bl_x0 = Vector3d(0, 0, half_pcd);
	Vector3d bl_x1 = Vector3d(0, -0.5 * half_pcd, 0.866025 * half_pcd);
	Vector3d bl_x2 = Vector3d(0, -0.866025 * half_pcd, 0.5 * half_pcd);
	Vector3d bl_x3 = Vector3d(0, -half_pcd, 0);
	Vector3d bl_x4 = Vector3d(0, -0.866025 * half_pcd, -0.5 * half_pcd);
	Vector3d bl_x5 = Vector3d(0, -0.5 * half_pcd, -0.866025 * half_pcd);
	Vector3d bl_x6 = Vector3d(0, 0, -half_pcd);

	EXPECT_NEAR((this->b4p.BL[0].x - bl_x0).norm(), 0, bl_x0.norm() *err);
	EXPECT_NEAR((this->b4p.BL[1].x - bl_x1).norm(), 0, bl_x1.norm() *err);
	EXPECT_NEAR((this->b4p.BL[2].x - bl_x2).norm(), 0, bl_x2.norm() *err);
	EXPECT_NEAR((this->b4p.BL[3].x - bl_x3).norm(), 0, bl_x3.norm() *err);
	EXPECT_NEAR((this->b4p.BL[4].x - bl_x4).norm(), 0, bl_x4.norm() *err);
	EXPECT_NEAR((this->b4p.BL[5].x - bl_x5).norm(), 0, bl_x5.norm() *err);
	EXPECT_NEAR((this->b4p.BL[6].x - bl_x6).norm(), 0, bl_x6.norm() *err);

	// 玉速度の検証
	double cos_alp = 0.86602540378443864676372317075294;
	double omega_c = (1 - FI.BL[0].dia * FI.cos_alp0 / FI.ballpcd * 0.5) * FI.omegair * 0.5; // 保持器公転数[min-1]
	double bl_v = omega_c * half_pcd; // 玉速度（スカラー）
	Vector3d bl_v0 = Vector3d(0, -bl_v, 0);
	Vector3d bl_v1 = Vector3d(0, -0.866025 * bl_v, -0.5 * bl_v);
	Vector3d bl_v2 = Vector3d(0, -0.5 * bl_v, -0.866025 * bl_v);
	Vector3d bl_v3 = Vector3d(0, 0, -bl_v);
	Vector3d bl_v4 = Vector3d(0, 0.5 * bl_v, -0.866025 * bl_v);
	Vector3d bl_v5 = Vector3d(0, 0.866025 * bl_v, -0.5 * bl_v);
	Vector3d bl_v6 = Vector3d(0, bl_v, 0);

	EXPECT_NEAR((this->b4p.BL[0].v - bl_v0).norm(), 0, bl_v0.norm() *err);
	EXPECT_NEAR((this->b4p.BL[1].v - bl_v1).norm(), 0, bl_v1.norm() *err);
	EXPECT_NEAR((this->b4p.BL[2].v - bl_v2).norm(), 0, bl_v2.norm() *err);
	EXPECT_NEAR((this->b4p.BL[3].v - bl_v3).norm(), 0, bl_v3.norm() *err);
	EXPECT_NEAR((this->b4p.BL[4].v - bl_v4).norm(), 0, bl_v4.norm() *err);
	EXPECT_NEAR((this->b4p.BL[5].v - bl_v5).norm(), 0, bl_v5.norm() *err);
	EXPECT_NEAR((this->b4p.BL[6].v - bl_v6).norm(), 0, bl_v6.norm() *err);

	// 角速度の検証
	double wx = -(FI.ballpcd / FI.BL[0].dia / 0.5 / cos_alp - FI.BL[0].dia * cos_alp * 0.5 / FI.ballpcd) * 0.5 * FI.omegair;
	Vector3d bl_w0 = Vector3d(wx, 0, 0);
	EXPECT_NEAR((this->b4p.BL[0].w - bl_w0).norm(), 0, bl_w0.norm() *err);
	EXPECT_NEAR((this->b4p.BL[1].w - bl_w0).norm(), 0, bl_w0.norm() *err);
	EXPECT_NEAR((this->b4p.BL[2].w - bl_w0).norm(), 0, bl_w0.norm() *err);
	EXPECT_NEAR((this->b4p.BL[3].w - bl_w0).norm(), 0, bl_w0.norm() *err);
	EXPECT_NEAR((this->b4p.BL[4].w - bl_w0).norm(), 0, bl_w0.norm() *err);
	EXPECT_NEAR((this->b4p.BL[5].w - bl_w0).norm(), 0, bl_w0.norm() *err);
	EXPECT_NEAR((this->b4p.BL[6].w - bl_w0).norm(), 0, bl_w0.norm() *err);

	//玉と内輪の座標・速度を取得
	this->b4p.save(FO);
	Vector3d BLx0(FO.BL[0].x);
	EXPECT_NEAR((BLx0 - bl_x0).norm(), 0, bl_x0.norm() *err);
	Vector3d BLv0(FO.BL[0].v);
	EXPECT_NEAR((BLv0 - bl_v0).norm(), 0, bl_v0.norm() *err);
	Vector3d BLw0(FO.BL[0].w);
	EXPECT_NEAR((BLw0 - bl_w0).norm(), 0, bl_w0.norm() *err);
	Vector3d IRx0(FO.IR.x);
	Vector3d ir_x = Vector3d(1, 0, 0);
	EXPECT_NEAR((IRx0 - ir_x).norm(), 0, ir_x.norm() *err);
	Vector3d IRv0(FO.IR.v);
	Vector3d ir_v = Vector3d(0, 0, 0);
	EXPECT_NEAR((IRv0 - ir_v).norm(), 0, ir_v.norm() *err);
	Vector3d IRw0(FO.IR.w);
	Vector3d ir_w = Vector3d(31.4 * 0.7071067, 0, 31.4 * 0.7071067);
	EXPECT_NEAR((IRw0 - ir_w).norm(), 0, ir_w.norm() *err);

	return;
};

// 変数yの入出力をテスト（単位時間やクォータニオンの計算が正しくできているか確認）
TEST_F(B4P_BearingTest, get_set_y) {
	int ballnum = 10;
	int nX = 13 * ballnum + 26;						//玉・内外輪・保持器各要素の状態量（座標・速度・加速度）
	double*y = new double[nX];
	this->test_init(ballnum);
	for (int i = 0; i < ballnum + 2; i++)
		for (int j = 0; j < 13; j++)
			y[i * 13 + j] = (i + 1) * 0.001;

	Rigid::l = 0.4;
	Rigid::t = 2.0;
	b4p.set_y(y);

	double err = 0.01;		// 許容誤差
	int i = 0;
	Vector3d bl_x = this->b4p.BL[i].x;
	Vector3d bl_v = this->b4p.BL[i].v;
	Quaterniond bl_q = this->b4p.BL[i].q;
	Vector3d bl_w = this->b4p.BL[i].w;
	Vector3d bl_x_ex = Vector3d::Ones() * (i + 3) * 0.001 * 0.4;
	Vector3d bl_v_ex = Vector3d::Ones() * (i + 3) * 0.001 * 0.4 * 0.5;
	Quaterniond bl_q_ex = Quaterniond(0.5, 0.5, 0.5, 0.5);
	Vector3d bl_w_ex = Vector3d::Ones() * (i + 3) * 0.001 * 0.5;
	EXPECT_NEAR((bl_x - bl_x_ex).norm(), 0, bl_x_ex.norm()*err);
	EXPECT_NEAR((bl_v - bl_v_ex).norm(), 0, bl_v_ex.norm()*err);
	double diff = Numeric::Square(bl_q.x() - bl_q_ex.x()) + Numeric::Square(bl_q.y() - bl_q_ex.y())
		+ Numeric::Square(bl_q.z() - bl_q_ex.z()) + Numeric::Square(bl_q.w() - bl_q_ex.w());
	EXPECT_NEAR(diff, 0, diff*err);
	EXPECT_NEAR((bl_w - bl_w_ex).norm(), 0, bl_w_ex.norm()*err);

	b4p.get_y(y);
	for (int j = 0; j < 13; j++) {
		double y_ex = (i + 1) * 0.001;
		if (j >= 6 && j <= 9) {
			y_ex = 0.5;
		}
		EXPECT_NEAR(y[i * 13 + j], y_ex, y_ex*err);
	}
	return;
}
//
//// 純アキシアル荷重条件下（接触角30度）で玉自転軸を100通り計算し，自転軸角30度で最もトルクが低くなることを確認するテスト．
//TEST_F(B4P_BearingTest, get_dydt1) {
//	int ballnum = 1;
//	this->test_init(ballnum);
//
//	double dx0 = 1e-5;
//	double x0[3], ax0[3];
//	x0[0] = dx0;  
//	x0[1] = 0.0;  
//	x0[2] = 0.0;
//	ax0[0] = 1.0; 
//	ax0[1] = 0.0; 
//	ax0[2] = 0.0;
//	this->FI.omegair = 31.4;
//	this->FI.omegaor = 0.0;
//
//	B4P_Bearing bearing;
//	bearing.init(this->FI, true, true, x0, ax0);
//	Rigid::g = Vector3d::Zero();
//	Rigid::l = 1.0;
//	Rigid::t = 1.0;
//	bearing.F_load = Vector3d::Zero();
//	bearing.T_load = Vector3d::Zero();
//
//	bearing.BL[0].x = Vector3d(dx0 / 2, 0, this->FI.ballpcd / 2);
//
//	int nXstf = 2 * ballnum + 5;		// 玉・内外輪・保持器各要素の状態量（座標・速度・加速度）
//	double*Xstf0 = new double[nXstf];
//	double*Fstf0 = new double[nXstf];
//	bearing.get_Xstf(Xstf0);
//	bearing.get_Fstf(Fstf0);
//
//	return;
//}
