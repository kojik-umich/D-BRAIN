#include "pch.h"
#include "BS_FileIn_stab.h"

class BS_CylinderTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	BS_FileIn_stab IN;
	BS_Shaft ST;

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
			this->IN.shaft[0].spiral[0].groove[i].r = 0.2;
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

		Rigid::l = 1;
		Rigid::t = 1;

		this->ST.init_pos(v0, w0);

		this->ST.x = Vector3d(0.0, 0.0, 0.0);
		this->ST.set_dx();
	}


	virtual void TearDown() {
	}
};

// 座標を螺旋座標系に変換して戻す（詳細なテストはSpiralクラスで実施済みであるため，単純なもののみ実施）
TEST_F(BS_CylinderTest, to_etacoord0) {
	this->init();

	double th = 1e1;

	Vector3d eta0(th, 1e-3, 2e-3);
	Vector3d x0 = this->ST.CY->to_inertialcoord(0, eta0);

	Vector3d eta1 = this->ST.CY->to_etacoord(0, x0);

	double delta = (eta0 - eta1).norm();

	EXPECT_GT(1e-2, delta);
};

// 速度ベクトルを螺旋座標系に変換して戻す
TEST_F(BS_CylinderTest, to_etavelocity0) {
	this->init();

	double th = Numeric::pi / 3;
	Matrix3d xyz2eta = this->ST.CY->get_xyz2eta(0, th);
	Vector3d v(1, 1e-3, 2e-3);
	Vector3d v_eta = this->ST.CY->to_etavector(v, xyz2eta);

	Vector3d v_ine = this->ST.CY->to_inertialvelocity(v_eta, xyz2eta);

	double delta = (v - v_ine).norm();

	EXPECT_GT(1e-2, delta);
};

// 位相角0°の位置で玉と溝が接触した際の接触位置を導出
TEST_F(BS_CylinderTest, calc_slice_test1) {
	this->init();

	int is = 0, ig = 0; // 条番号，溝番号
	double th = 0;	//	螺旋位相角
	Vector3d p(0, 0.8 - 1e-6, 0);	// 接触点[m]
	Vector3d t(0.8, 0, 2 * Numeric::pi * 1.0);	// 螺旋進行ベクトル[m]
	Vector3d xai = t / t.norm();				// 螺旋進行方向単位ベクトル[-]
	double a = 1e-4;							// 接触楕円長径（適当）
	int ms = 21;								// スライス分割数							
	double bl_r = 0.15;							// 玉半径[m](適当)
	Vector3d ps[21];

	this->ST.CY->calc_slice(is, ig, th, p, xai, a, ms, bl_r, ps);

	// (1) 接触楕円の端から端と，トーラス中心から接触点へのベクトルが直行しているかのかの確認
	Vector3d op = Vector3d(0, 1, 0);
	double ip1 = op.normalized().dot((ps[20] - ps[0]).normalized());
	EXPECT_NEAR(ip1, 0, 1e-6);

	// (2) 接触楕円の端から端と，螺旋進行方向が直行しているかの確認
	double ip2 = xai.dot((ps[20] - ps[0]).normalized());
	EXPECT_NEAR(ip2, 0, 1e-6);

	// (3) 接触楕円の端から端の距離が，接触楕円の長さに等しいか確認
	// (スライス一つ分だけ短くなるため，加味して計算)
	double ip3 = (ps[20] - ps[0]).norm();
	EXPECT_NEAR(ip3, a * 2 / 21.0 * 20.0 , 1e-6);

	// (4) 接触楕円の片端から溝R中心への距離が両方等しいかの確認
	Vector3d Og(0, 1.0, 0);	// 溝R中心[m]
	double l1 = (ps[20] - Og).norm();
	double l2 = (ps[0] - Og).norm();
	EXPECT_NEAR(l1, l2, 1e-6);

	// (5) 接触楕円の片端から溝R中心への距離より，接触点から中心への距離の方が大きいかの確認
	double l0 = (ps[20] - Og).norm();
	EXPECT_GE(l0, l1);

	// (6) 接触楕円の片端から溝R中心への距離が，ザックリ溝半径と等しいかの確認
	EXPECT_NEAR(l1, 0.2, 1e-5);

	return;
}
