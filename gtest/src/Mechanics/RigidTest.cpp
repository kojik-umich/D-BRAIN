#include "pch.h"

class RigidTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	Ball BL;

	virtual void SetUp() {

		Rigid::l = 1e-2;
		Rigid::t = 1e-3;
		Rigid::g = Vector3d(0.0, 0.0, 9.8);
	}
};

// 軸方向の指定メソッドのテスト．軸を指定したのち軸の方向を取得し，もとの値と一致しているかを確認．2通りの出力で計算．
TEST_F(RigidTest, setget_ax) {

	Vector3d ax0 = Vector3d(0.8, -0.6, 0.0);

	this->BL.set_ax(ax0);
	Vector3d return0_ = this->BL.get_ax();
	double   delta0_  = (return0_ - ax0).norm();

	EXPECT_GT(1e-2, delta0_);

	Vector3d ax1 = Vector3d(0.0, 0.8, 0.6);

	this->BL.set_ax(ax1);
	Vector3d return1_ = this->BL.get_ax();
	double   delta1_  = (return1_ - ax1).norm();

	EXPECT_GT(1e-2, delta1_);

	return;
};

// SI単位系のまま入力して出力する値の一致を確認するテスト．
TEST_F(RigidTest, setget_param0) {

	Vector3d x0 = Vector3d(0.0, 0.8, 0.6);
	Vector3d v0 = Vector3d(1.0, 2.8, 5.6);
	Quaterniond q0 = Quaterniond(0.5, -0.5, 0.5, -0.5);
	Vector3d w0 = Vector3d(0.7, 6.8, 2.6);
	Vector3d x1, v1, w1;
	Quaterniond q1;

	this->BL.set_param(x0, v0, q0, w0);
	this->BL.get_param(x1, v1, q1, w1);
	double   delta_x = (x1 - x0).norm();
	double   delta_v = (v1 - v0).norm();
	double   delta_q = abs(q1.w() - q0.w()) + abs(q1.x() - q0.x()) + abs(q1.y() - q0.y()) + abs(q1.z() - q0.z());
	double   delta_w = (w1 - w0).norm();

	EXPECT_GT(1e-2, delta_x);
	EXPECT_GT(1e-2, delta_v);
	EXPECT_GT(1e-2, delta_q);
	EXPECT_GT(1e-2, delta_w);

	return;
};

// 姿勢の時間微小変化量が正しく計算できているかを確認するテスト．変化前と後の軸方向をそれぞれ取り出し，角速度の定義通りに回転しているかを確認する．
TEST_F(RigidTest, get_dqdt0) {

	Quaterniond q0 = Quaterniond(1.0, 0.0, 0.0, 0.0);
	this->BL.q = q0;

	Vector3d w0 = Vector3d(0.0, 1.0, 0.0);
	this->BL.w = w0;

	Quaterniond dqdt = this->BL.get_dqdt();
	double dt0 = 1.0e-6;

	this->BL.q = Quaterniond(
		q0.w() + dqdt.w() * dt0,
		q0.x() + dqdt.x() * dt0,
		q0.y() + dqdt.y() * dt0,
		q0.z() + dqdt.z() * dt0
	);

	double th  = w0.y() * dt0;

	Vector3d ax0 = this->BL.get_ax();
	Vector3d ax1 = Vector3d(cos(th), 0.0, -sin(th));

	double   delta_x = (ax1 - ax0).norm();

	EXPECT_NEAR(delta_x, 0, 1e-6);

	return;
};

// 角速度の時間微小変化量が正しく計算できているかを確認するテスト．世間一般の有名な慣性モーメントの式から手計算し，値の一致を確認する．
TEST_F(RigidTest, get_dwdt0) {

	Vector3d w0 = Vector3d(0.0, 0.0, 0.0);
	this->BL.w = w0;

	Vector3d I0 = Vector3d(2.0, 3.0, 4.0);
	this->BL.I = I0;
	this->BL.I_inv = I0.cwiseInverse();

	Vector3d T0 = Vector3d(1.0, 0.0, 0.0);

	Vector3d dwdt = this->BL.get_dwdt(T0);
	double dt0 = 1.0;

	Vector3d dw0 = dwdt * dt0;
	Vector3d dw1 = Vector3d(0.5, 0.0, 0.0);

	double   delta_x = (dw1 - dw0).norm();

	EXPECT_GT(1e-2, delta_x);

	return;
};

// 物体の表面速度の計算が正しくできているかを確認するテスト．手計算との一致を確認．
TEST_F(RigidTest, surface_velocity0) {

	Vector3d x0 = Vector3d(1.0, 0.0, 0.0);
	this->BL.x = x0;

	Vector3d v0 = Vector3d(0.0, 2.0, 0.0);
	this->BL.v = v0;

	Vector3d w0 = Vector3d(0.0, 0.0, 5.0);
	this->BL.w = w0;

	Vector3d x1 = Vector3d(0.0, 0.0, 0.0);

	Vector3d V0 = this->BL.surface_velocity(x1);
	Vector3d V1 = Vector3d(0.0, -3.0, 0.0);

	double   delta_x = (V1 - V0).norm();

	EXPECT_GT(1e-2, delta_x);

	return;
};

// (1) 指定した単位長さ・単位時間で規格化された時間微分値の取り出しテスト．手計算との一致を確認．
TEST_F(RigidTest, get_dydt_test1) {

	Vector3d x0 = Vector3d(1.0, 0.0, 0.0);
	this->BL.x = x0;

	Vector3d v0 = Vector3d(1.0, 2.0, 3.0);
	this->BL.v = v0;

	Quaterniond q0 = Quaterniond::Identity();
	this->BL.q = q0;

	Vector3d w0 = Vector3d(4.0, 5.0, 6.0);
	this->BL.w = w0;

	Vector3d F0 = Vector3d(7.0, 8.0, 9.0);
	Vector3d T0 = Vector3d(10.0, 11.0, 12.0);

	this->BL.m = 1.0;
	this->BL.m_inv = 1.0;
	this->BL.I = Vector3d::Ones();
	this->BL.I_inv = Vector3d::Ones();

	double dydt0[13];
	this->BL.get_dydt(F0, T0, dydt0);

	Vector3d x1 = Vector3d(100.0, 0.0, 0.0);
	this->BL.x = x1;
	this->BL.get_dydt(F0, T0, dydt0);

	EXPECT_GT(1e-2, abs(Rigid::l / Rigid::t * dydt0[0] - 1.0));
	EXPECT_GT(1e-2, abs(Rigid::l / Rigid::t * dydt0[1] - 2.0));
	EXPECT_GT(1e-2, abs(Rigid::l / Rigid::t * dydt0[2] - 3.0));
	EXPECT_GT(1e-2, abs(Rigid::l / Rigid::t / Rigid::t * dydt0[3] - 7.0));
	EXPECT_GT(1e-2, abs(Rigid::l / Rigid::t / Rigid::t * dydt0[4] - 8.0));
	EXPECT_GT(1e-2, abs(Rigid::l / Rigid::t / Rigid::t * dydt0[5] - 9.0));
	EXPECT_GT(1e-2, abs(1 / Rigid::t * dydt0[6] - 0.0));
	EXPECT_GT(1e-2, abs(1 / Rigid::t * dydt0[7] - 2.0));
	EXPECT_GT(1e-2, abs(1 / Rigid::t * dydt0[8] - 2.5));
	EXPECT_GT(1e-2, abs(1 / Rigid::t * dydt0[9] - 3.0));
	EXPECT_GT(1e-2, abs(1 / Rigid::t / Rigid::t * dydt0[10] - 10.0));
	EXPECT_GT(1e-2, abs(1 / Rigid::t / Rigid::t * dydt0[11] - 11.0));
	EXPECT_GT(1e-2, abs(1 / Rigid::t / Rigid::t * dydt0[12] - 12.0));

	return;
};

// (2) 拘束条件が設定されている場合の挙動の確認
TEST_F(RigidTest, get_dydt_test2) {
	Vector3d x0 = Vector3d(1.0, 0.0, 0.0);
	this->BL.x = x0;
	Vector3d v0 = Vector3d(1.0, 2.0, 3.0);
	this->BL.v = v0;
	Quaterniond q0 = Quaterniond::Identity();
	this->BL.q = q0;
	Vector3d w0 = Vector3d(4.0, 5.0, 6.0);
	this->BL.w = w0;

	Vector3d F0 = Vector3d(7.0, 8.0, 9.0);
	Vector3d T0 = Vector3d(10.0, 11.0, 12.0);
	this->BL.set_mI(1.0, Vector3d::Ones());

	bool v_const[3], w_const[3];
	for (int i = 0; i < 3; i++) {
		v_const[i] = true;
		w_const[i] = true;
	}

	this->BL.set_const(v_const, w_const);
	double dydt0[13];
	this->BL.get_dydt(F0, T0, dydt0);

	EXPECT_NEAR(dydt0[3], 0, 1e-6);
	EXPECT_NEAR(dydt0[4], 0, 1e-6);
	EXPECT_NEAR(dydt0[5], 0, 1e-6);
	EXPECT_NEAR(dydt0[10], 0, 1e-6);
	EXPECT_NEAR(dydt0[11], 0, 1e-6);
	EXPECT_NEAR(dydt0[12], 0, 1e-6);
	return;
}

// クォータニオンによるベクトル回転の演算が正しく計算できているかのテスト．[x,y,z]ベクトルを軸に，2θだけ回転するクォータニオンが[cosθ, xsinθ, ysinθ, zsinθ]であることを手計算との一致により確認．
TEST_F(RigidTest, to_myvelocity0) {

	this->BL.v = Vector3d(2.0, 0.0, 0.0);

	// z軸正回り 90度回転姿勢．
	this->BL.q = Quaterniond(sqrt(0.5), 0.0, 0.0, sqrt(0.5));

	Vector3d a(0.0, 2.0, 0.0);
	Vector3d b = this->BL.to_myvector(a);
	Vector3d b_(2.0, 0.0, 0.0);

	double   db = (b - b_).norm();
	EXPECT_GT(1e-2, db);

	Vector3d a_ = this->BL.to_inevector(b);
	double   da = (a - a_).norm();
	EXPECT_GT(1e-2, da);

	Vector3d v0 = Vector3d(3.0, 0.0, 0.0);

	Vector3d V0 = this->BL.to_myvelocity(v0);
	Vector3d V1 = Vector3d(0.0, -1.0, 0.0);

	double   dv0 = (V1 - V0).norm();
	EXPECT_GT(1e-2, dv0);

	Vector3d v0_ = this->BL.to_inevelocity(V0);

	double   dv1 = (v0_ - v0).norm();
	EXPECT_GT(1e-2, dv1);

	return;
};

// 引き続きクォータニオンの確認および，部材のy軸，z軸が正しく計算できているかの確認．手計算との一致を確認．
TEST_F(RigidTest, get_ay_az) {

	// x軸正回り 90度回転姿勢．
	this->BL.q = Quaterniond(sqrt(0.5), sqrt(0.5), 0.0, 0.0);

	Vector3d a = this->BL.get_ay();
	Vector3d a_(0.0, 0.0, 1.0);

	double   da = (a - a_).norm();
	EXPECT_GT(1e-2, da);

	Vector3d b = this->BL.get_az();
	Vector3d b_(0.0, -1.0, 0.0);

	double   db = (b - b_).norm();
	EXPECT_GT(1e-2, db);

	Vector3d T(1.0, 2.0, 3.0);
	Vector3d c = this->BL.remove_ax(T);
	Vector3d c_(0.0, 2.0, 3.0);
	double   dc = (c - c_).norm();
	EXPECT_GT(1e-2, dc);

	return;
};

// 中心位置が原点にない部材に外力をかけた際のトルクの計算が手計算と合っているかの確認．2通りの計算で実施．
TEST_F(RigidTest, calc_Torque0) {

	this->BL.x = Vector3d(2.0, 0.0, 0.0);
	Vector3d p(4.0, 0.0, 0.0);
	Vector3d F(0.0, 3.0, 0.0);

	Vector3d T = this->BL.calc_Torque(p, F);
	Vector3d T_(0.0, 0.0, 6.0);

	double   dT = (T - T_).norm();
	EXPECT_GT(1e-2, dT);

	Vector3d u(1.0, 3.0, 0.0);
	Vector3d eT = this->BL.calc_TorqueDirection(p, u);
	Vector3d eT_(0.0, 0.0, 1.0);

	double   deT = (eT - eT_).norm();
	EXPECT_GT(1e-2, deT);

	return;
};

// 重力ベクトルを正しく取得できるか確認のテスト．
TEST_F(RigidTest, get_mg0) {

	this->BL.m = 10;

	Vector3d mg = this->BL.get_mg();
	Vector3d mg_(0.0, 0.0, 98);

	double dmg = (mg - mg_).norm();
	EXPECT_GT(1e-2, dmg);

	return;
};

// 出力用の変数が正しく計算できるかのテスト．
TEST_F(RigidTest, make_OutParam0) {

	Vector3d x, v, w, F, T;
	x = v = w = F = T = Vector3d::Ones();
	Quaterniond q;
	q.setIdentity();
	this->BL.set_param(x, v, q, w);
	this->BL.set_FT(F, T);


	double *_x, *_v, *_q, *_w, *_ax, *_F, *_T;
	_x = new double[3];
	_v = new double[3];
	_w = new double[3];
	_q = new double[4];
	_ax = new double[3];
	_F = new double[3];
	_T = new double[3];
	this->BL.save(_x, _v, _q, _w, _ax, _F, _T);

	for (int i = 0; i < 3; i++) {
		EXPECT_NEAR(_x[i], 1, 1e-6);
		EXPECT_NEAR(_v[i], 1, 1e-6);
		EXPECT_NEAR(_w[i], 1, 1e-6);
		EXPECT_NEAR(_F[i], 1, 1e-6);
		EXPECT_NEAR(_T[i], 1, 1e-6);
	}

	EXPECT_NEAR(_q[0], 0, 1e-6); EXPECT_NEAR(_q[1], 0, 1e-6); 
	EXPECT_NEAR(_q[2], 0, 1e-6); EXPECT_NEAR(_q[3], 1, 1e-6);
	EXPECT_NEAR(_ax[0], 1, 1e-6); EXPECT_NEAR(_ax[1], 0, 1e-6);
	EXPECT_NEAR(_ax[2], 0, 1e-6);
	return;
};

// 残りの，計算がほとんど要らないけど一応確認しておくメソッド群．
TEST_F(RigidTest, all_rest) {

	Vector3d F(0.0, 0.0, 0.0);
	Vector3d T(0.0, 0.0, 0.0);

	this->BL.set_FT(F, T);

	Vector3d x, v, w;
	x = v = w = Vector3d::Ones();
	Quaterniond q;
	q.setIdentity();
	this->BL.set_y(x, v, q, w);

	Vector3d x_, v_, w_;
	Quaterniond q_;
	this->BL.get_y(x_, v_, q_, w_);

	double dy = (x - x_).norm() + (v - v_).norm() + (w - w_).norm();
	EXPECT_GT(1e-2, dy);

	return;
};

// 慣性座標系⇔物体座標系の変換，ただし物体の軸方向が(√3/2, 1/2, 0)で，座標が(1e-2, 2e-2, 3e-2)のときを考える．
TEST_F(RigidTest, to_mycoord_test1) {
	
	Vector3d ax(0.866025, 0.5, 0);
	this->BL.set_ax(ax);								// set_ax は検証済として扱う．
	this->BL.x = Vector3d(1e-2, 2e-2, 3e-2);
	
	// 慣性座標系から物体座標系への変換
	Vector3d x_in(2e-2, 9e-2, 0);						// 入力
	Vector3d x_ex1(4.366025e-2, 5.5621778e-2, -3e-2);	// 答え(予測値)
	Vector3d x_out1 = this->BL.to_mycoord(x_in);		// 答え(単体テスト)
	EXPECT_NEAR((x_out1 - x_ex1).norm(), 0, 1e-6);
	
	// 物体座標系から慣性座標系への変換
	// （最初の入力に戻るか判定）
	Vector3d x_out2 = this->BL.to_inecoord(x_out1);		// 答え(単体テスト)
	EXPECT_NEAR((x_out2 - x_in).norm(), 0, 1e-6);

	return;
}

