#include "pch.h"
#include "B4P_StabIn.h"

class B4P_OuterRingTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	B4P_OuterRing  OR;
	B4P_InnerRing  IR;
	double D;
	virtual void SetUp() {
		B4P_StabIn FI;

		// 内輪物性値
		this->D = 0.00635;				// 玉径[m]
		FI.ballpcd = 0.0355;
		double clr = 0;
		FI.IR.Rox[0] = 0.0000725;
		FI.IR.Rox[1] = -0.0000725;
		FI.IR.R[0] = 0.00332105;
		FI.IR.R[1] = 0.00332105;
		// 溝R_pcd (= 0.035753569339629220)
		FI.IR.Rod[0] = FI.ballpcd + (FI.IR.R[0] - D * 0.5) * 0.866025403784439;
		FI.IR.Rod[1] = FI.ballpcd + (FI.IR.R[1] - D * 0.5) * 0.866025403784439;
		FI.IR.hedge[0];
		FI.IR.hedge[1];
		FI.IR.rms = 0.00000006;
		FI.IR.m = 1;
		FI.IR.E = 207760000000;
		FI.IR.por = 0.29;
		FI.IR.Ix, FI.IR.Iyz, FI.IR.Iyz;

		// 潤滑理論
		FI.TB.rollingresistance;
		this->IR.init(FI.IR, FI.ballpcd);

		// 外輪物性値
		FI.OR.Rox[0] = 0.0000725;
		FI.OR.Rox[1] = -0.0000725;
		FI.OR.R[0] = 0.00332105;
		FI.OR.R[1] = 0.00332105;
		// 溝R_pcd (= 0.035246430660370774)
		FI.OR.Rod[0] = FI.ballpcd - (FI.OR.R[0] - D * 0.5) * 0.866025403784439;
		FI.OR.Rod[1] = FI.ballpcd - (FI.OR.R[1] - D * 0.5) * 0.866025403784439;
		FI.OR.hedge[0];
		FI.OR.hedge[1];
		FI.OR.rms = 0.00000006;
		FI.OR.m = 1;
		FI.OR.E = 207760000000;
		FI.OR.por = 0.29;
		FI.OR.Ix, FI.OR.Iyz, FI.OR.Iyz;

		this->OR.init(FI.OR, FI.ballpcd);

		Vector3d x = Vector3d(0, 0, 0);
		Vector3d v = Vector3d(0, 0, 0);
		Quaterniond q = Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）
		Vector3d w = Vector3d(0, 0, 0);

		this->OR.set_y(x, v, q, w);
		this->IR.set_y(x, v, q, w);
	}

	virtual void TearDown() {

	}
};

TEST_F(B4P_OuterRingTest, Nothing) {

	this->OR.Nothing();
	EXPECT_GT(1e-2, 0.0);
};


// (1) 各スライス片の中心と溝R中心や各種方向ベクトルの相対関係が正しいか判定
TEST_F(B4P_OuterRingTest, calc_slice_test1) {

	Vector3d p(this->OR.GV[0].Rx, 0, this->OR.GV[0].r + this->OR.GV[0].Rr);
	double a = 1e-4;
	int ms = 21;
	Vector3d ps[21];

	this->OR.calc_slice(p, a, ms, 0, this->D, ps);
	EXPECT_NEAR((ps[20] - ps[0]).norm(), 2 * a, 1e-1 * a) << "接触楕円の端から端の距離が，接触楕円の長さに等しいか確認";

	// 接触楕円の端から端と，溝R中心から接触点へのベクトルが直行しているかのかの確認
	Vector3d OP(0, 0, this->OR.GV[0].r);
	double ip1 = OP.dot((ps[20] - ps[0]));
	EXPECT_NEAR(ip1, 0 , 1e-6);

	// 接触楕円の端から端と，周方向ベクトルが直交しているかの確認
	Vector3d xai(0, 1, 0);
	double ip2 = xai.dot(ps[20] - ps[0]);
	EXPECT_NEAR(ip2, 0, 1e-6);

	// スライス片中心から溝R中心への距離がそれぞれ等しいかの確認
	Vector3d O(this->OR.GV[0].Rx, 0, this->OR.GV[0].Rr);
	double l1 = (ps[20] - O).norm();
	double l2 = (ps[0] - O).norm();
	EXPECT_NEAR(l1, l2, 1e-6);

	//  接触楕円の片端から溝R中心への距離 l1より接触点から溝R中心への距離 l0 が大きいかの確認
	// (仕様検討のため保留)
	double l0 = (ps[10] - O).norm();
	EXPECT_GE(l0, l1);

	// 接触楕円の片端から溝中心への距離が，溝半径と1%以内の差にとどまっているかの確認
	Vector3d Gc(this->OR.GV[0].Rx, 0, this->OR.GV[0].Rr);
	double OP0 = (ps[20] - Gc).norm();
	double r = this->OR.GV[0].r;
	EXPECT_NEAR(l1, r, 1e-5) << "接触楕円の片端から溝中心への距離が，溝半径と1%以内の差にとどまっているかの確認";
	return;
};


// (2) 二次関数近似ではなく円弧の式を用いて計算，比較
TEST_F(B4P_OuterRingTest, calc_slice_test2) {
	// 接触角30°の時の接触点
	Vector3d Ro(this->OR.GV[0].Rx, 0, this->OR.GV[0].Rr);	// 着目する溝R中心
	double theta = Numeric::pi * 0.166666;					// 接触角
	Vector3d p = Ro + this->OR.GV[0].r * Vector3d(sin(theta), 0, cos(theta));

	// プログラムの式で計算
	double a = 1e-6;					// 接触楕円[m]
	int ms = 21;						// 分割数
	Vector3d *ps = new Vector3d[ms];	// 接触点座標[m]
	this->OR.calc_slice(p, a, ms, 0, this->D, ps);


	// 円の式を用いて計算
	Vector3d *ps_ex = new Vector3d[ms];	// 接触点座標（期待値）[m]
	for (int j = 0; j < ms; j++) {
		double z = -((ms - 1) * 0.5 - j) * (2 * a / ms * sin(theta)) + p.z();
		double x = sqrt(this->OR.GV[0].r * this->OR.GV[0].r - (z - Ro.z()) * (z - Ro.z())) + Ro.x();
		ps_ex[j] = Vector3d(x, 0, z);
	}

	for (int j = 0; j < ms; j++) {
		double err = 1e-8;
		EXPECT_NEAR(ps[j].x(), ps_ex[j].x(), err);
		EXPECT_NEAR(ps[j].y(), ps_ex[j].y(), err);
		EXPECT_NEAR(ps[j].z(), ps_ex[j].z(), err);
	}
	delete[] ps;
	return;
}

// (3) 内輪と外輪で玉を接触角30°になるように挟んだ時，内輪と外輪のスライス片が並行になるか検証．
TEST_F(B4P_OuterRingTest, calc_slice_test3) {
	// 内輪x座標 = 2 * (溝R半径 - 玉半径) * sin30° - Rx * 2
	double bl_r = 0.00635 * 0.5;
	double Ox = 2 * (this->IR.GV[0].r - bl_r) * 0.5 - (this->IR.GV[0].Rx) * 2;
	Vector3d x_out	= Vector3d(0, 0, 0);
	Vector3d x_in	= Vector3d(Ox, 0, 0);
	Vector3d v		= Vector3d(0, 0, 0);
	Quaterniond q	= Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）
	Vector3d w		= Vector3d(0, 0, 0);

	this->OR.set_y(x_out, v, q, w);
	this->IR.set_y(x_in, v, q, w);

	double a = 1e-6;								// 接触楕円[m]
	int ms = 21;									// 分割数
	Vector3d *ps_out = new Vector3d[ms];			// 接触点座標[m]
	Vector3d p_out = Vector3d(this->OR.GV[1].Rx, 0, this->OR.GV[1].Rr) + this->OR.GV[1].r * Vector3d(0.5, 0, 0.866025403784439);
	this->OR.calc_slice(p_out, a, ms, 1, this->D, ps_out);
	
	Vector3d p_in_other = Vector3d(Ox + this->IR.GV[0].Rx, 0, this->IR.GV[0].Rr) + this->IR.GV[0].r * Vector3d(-0.5, 0, -0.866025403784439);
	Vector3d *ps_in = new Vector3d[ms];				// 接触点座標[m]
	Vector3d p_in  = p_out + 2 * bl_r * Vector3d(-0.5, 0, -0.866025403784439);
	this->IR.calc_slice(p_in_other, a, ms, 0, this->D, ps_in);
	Vector3d pdist_out = ps_out[20] - ps_out[0];
	Vector3d pdist_in = ps_in[20] - ps_in[0];

	Vector3d cross = pdist_out.normalized().cross(pdist_in.normalized());
	EXPECT_NEAR(cross.norm(), 0, 1e-3);
	delete[] ps_out, ps_in;

	return;
}

// 慣性座標系から溝直交断面座標系への変換
TEST_F(B4P_OuterRingTest, ine_to_XZvector_1) {

	double sqrt3 = 1.7320508075688772935274463415059;
	Vector3d a(1e-3, 10e-3, sqrt3 * 1e-2);
	Vector3d thXZ = this->OR.ine_to_XZcoord(a);

	Vector3d thXZ_ans(-Numeric::pi * 0.16666666, 1e-3, 2.25e-3);
	EXPECT_NEAR((thXZ - thXZ_ans).norm(), 0, 1e-6);
	Matrix3d xyz2XYZ = this->OR.get_xyz2XYZ(thXZ[0]);

	Vector3d b(0, sqrt3, -1.0);
	Vector3d Gx = this->OR.ine_to_XZvector(b, xyz2XYZ);
	Vector3d Gx_ans(0, 2.0, 0);
	EXPECT_NEAR((Gx - Gx_ans).norm(), 0, 1e-6);

	return;
}


class B4P_InnerRingTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	B4P_InnerRing  IR;

	virtual void SetUp() {
	}

	virtual void TearDown() {

	}
};

TEST_F(B4P_InnerRingTest, Nothing) {

	this->IR.Nothing();
	EXPECT_GT(1e-2, 0.0);
};

// (1) 指定した単位長さ・単位時間で規格化された時間微分値の取り出しテスト．手計算との一致を確認．
TEST_F(B4P_InnerRingTest, get_dydt__1) {
	for (int i = 0; i < 3; i++) {
		IR.x_const[i] = false;
		IR.Rx_const[i] = false;
	}

	Vector3d x0 = Vector3d(1.0, 0.0, 0.0);
	this->IR.x = x0;

	Vector3d v0 = Vector3d(1.0, 2.0, 3.0);
	this->IR.v = v0;

	Quaterniond q0 = Quaterniond::Identity();
	this->IR.q = q0;

	Vector3d w0 = Vector3d(4.0, 5.0, 6.0);
	this->IR.w = w0;

	Vector3d F0 = Vector3d(7.0, 8.0, 9.0);
	Vector3d T0 = Vector3d(10.0, 11.0, 12.0);

	this->IR.m = 1.0;
	this->IR.m_inv = 1.0;
	this->IR.I = Vector3d::Ones();
	this->IR.I_inv = Vector3d::Ones();

	double dydt0[13];
	this->IR.get_dydt_(F0, T0, dydt0);

	Vector3d x1 = Vector3d(100.0, 0.0, 0.0);
	this->IR.x = x1;
	this->IR.get_dydt_(F0, T0, dydt0);

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
	EXPECT_GT(1e-2, abs(1 / Rigid::t / Rigid::t * dydt0[10] - 0.0));
	EXPECT_GT(1e-2, abs(1 / Rigid::t / Rigid::t * dydt0[11] - 11.0));
	EXPECT_GT(1e-2, abs(1 / Rigid::t / Rigid::t * dydt0[12] - 12.0));

	return;
};

// (2) 拘束条件が設定されている場合の挙動の確認
TEST_F(B4P_InnerRingTest, get_dydt__2) {
	for (int i = 0; i < 3; i++) {
		IR.x_const[i] = true;
		IR.Rx_const[i] = true;
	}
	Vector3d x0 = Vector3d(1.0, 0.0, 0.0);
	this->IR.x = x0;
	Vector3d v0 = Vector3d(1.0, 2.0, 3.0);
	this->IR.v = v0;
	Quaterniond q0 = Quaterniond::Identity();
	this->IR.q = q0;
	Vector3d w0 = Vector3d(4.0, 5.0, 6.0);
	this->IR.w = w0;

	Vector3d F0 = Vector3d(7.0, 8.0, 9.0);
	Vector3d T0 = Vector3d(10.0, 11.0, 12.0);
	this->IR.set_mI(1.0, Vector3d::Ones());

	double dydt0[13];
	this->IR.get_dydt_(F0, T0, dydt0);

	EXPECT_NEAR(dydt0[3], 0, 1e-6);
	EXPECT_NEAR(dydt0[4], 0, 1e-6);
	EXPECT_NEAR(dydt0[5], 0, 1e-6);
	EXPECT_NEAR(dydt0[10], 0, 1e-6);
	EXPECT_NEAR(dydt0[11], 0, 1e-6);
	EXPECT_NEAR(dydt0[12], 0, 1e-6);
	return;
}
