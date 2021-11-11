#include "pch.h"

class BallTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
	Ball  BL0, BL1;
	BallBallPair BBP;

	virtual void SetUp() {

		double dia = 1.0;
		double E   = 1.0;
		double por = 1.0;
		double den = 1.0;
		double rms = 1.0;
		bool v_const[3], w_const[3];
		for (int i = 0; i < 3; i++) {
			v_const[i] = false;
			w_const[i] = false;
		}

		BBP.link(&this->BL0, &this->BL1);

		this->BL0.init(dia, E, por, den, rms, v_const, w_const);
		this->BL1.init(dia, E, por, den, rms, v_const, w_const);

		this->BL0.Nothing();
	}
};

TEST_F(BallTest, calc_Contact_test1) {

	Vector3d x0 = Vector3d(1.0, 0.0, 0.0);
	Vector3d v0 = Vector3d(1.0, 1.0, 0.0);
	Quaterniond q0 = Quaterniond(1.0, 0.0, 0.0, 0.0);
	Vector3d w0 = Vector3d(1.0, 0.0, 0.0);
	this->BL0.set_param(x0, v0, q0, w0);

	double dx, dv;
	Vector3d edir;
	// 0:玉中心と接触点が完全に重なるケース．例外処理に飛ぶ．
	Vector3d Zero = Vector3d::Zero();
	bool return0 = this->BL0.calc_Contact(x0, Zero, dx, dv, edir);
	EXPECT_FALSE(return0);

	// 1:玉中心と接触点の距離が玉半径より大きいケース．例外処理に飛ぶ．
	Vector3d x1 = Vector3d(1.6, 0.0, 0.0);
	bool return1 = this->BL0.calc_Contact(x1, Zero, dx, dv, edir);
	EXPECT_FALSE(return1);

	// 2:玉に接触点が食い込むケース．手計算との一致を確認．
	Vector3d x2 = Vector3d(1.4, 0.0, 0.0);
	Vector3d v2 = Vector3d(0.9, 1.0, 0.0);
	bool return2 = this->BL0.calc_Contact(x2, v2, dx, dv, edir);
	Vector3d edir_ans = Vector3d(1, 0.0, 0.0);
	EXPECT_TRUE(return2);
	EXPECT_NEAR(dx, 0.1, 1e-3);
	EXPECT_NEAR(dv, 0.1, 1e-3);
	EXPECT_NEAR((edir - edir_ans).norm(), 0, 1e-3);

	return;
};

// (1) 玉法線成分が削除できているか確認
TEST_F(BallTest, remove_Normal_test1) {
	Vector3d x0 = Vector3d(1.0, 0.0, 0.0);
	Vector3d v0 = Vector3d(1.0, 1.0, 0.0);
	Quaterniond q0 = Quaterniond(1.0, 0.0, 0.0, 0.0);
	Vector3d w0 = Vector3d(1.0, 0.0, 0.0);
	this->BL0.set_param(x0, v0, q0, w0);

	Vector3d p(0.5, 0, 0);
	Vector3d a(1.0, 4.0, 0.5);
	Vector3d return1 = this->BL0.remove_Normal(p, a);
	Vector3d answer1(0.0, 4.0, 0.5);
	EXPECT_NEAR((return1 - answer1).norm(), 0, 1e-3);

	return;
}
