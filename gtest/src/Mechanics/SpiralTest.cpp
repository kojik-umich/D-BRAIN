#include "pch.h"
/*******************************************************************************
!								"SpiralTest.cpp"
!												2021/4/2	[Core-T]	楠崎，柏原
!	直交座標系⇔螺旋座標系の変換をテストしている．
!	テスト項目は
!	(1) "座標"を直交座標系⇒螺旋座標系に変換
!	(2) 螺旋座標に変換した"座標"をもとの直交座標系に変換できるか
!	(3) "ベクトル"を直交座標系⇒螺旋座標系に変換
!	(4) (3)の変換行列が回転行列の性質を満たしているか判定
!	
!	※B:\ソースコード他\D-BRAIN\docs\26_D-BS\03_評価\01_単体テスト\SpiralTest.pptxに補足説明を記載
!*******************************************************************************/


class SpiralTest : public ::testing::Test {
public:
	Spiral SP;
	
	virtual void SetUp() {
	}

	// テストデータ1
	void init_data01() {
		double alp = Unit::deg2rad(60);// 初期位相角[rad](30°)
		double l = 0.6;
		double r = 0.8;

		this->SP.init(alp, l, r);
	}

	// テストデータ2
	void init_data02() {
		double alp = Unit::deg2rad(0);	// 初期位相角[rad](0°)
		double l = 0.06;
		double r = 0.1;

		this->SP.init(alp, l, r);
	}
};
// (1)-1 "座標"を直交座標系⇒螺旋座標系に変換，手計算と比較
TEST_F(SpiralTest, x_to_eta_01) {
	this->init_data01();
	// 点A(x=0から300°分移動した位置)
	Vector3d x_A = Vector3d(0.5, 0.8 , 0);
	Vector3d eta_A = this->SP.to_eta2(x_A);
	Vector3d eta_A_ex = Vector3d(Unit::deg2rad(300), 0, 0);
	EXPECT_NEAR(eta_A.x(), eta_A_ex.x(), 1e-3);
	EXPECT_NEAR(eta_A.y(), eta_A_ex.y(), 1e-7);
	EXPECT_NEAR(eta_A.z(), eta_A_ex.z(), 1e-7);
	Vector3d x1 = this->SP.to_xyz(eta_A);
	Vector3d x2 = this->SP.to_xyz(eta_A_ex);
	Vector3d zeta = Vector3d(0.992951092, 0, -0.118524806); // ζの単位ベクトル

	// 点B(ζ方向に1mm移動)
	Vector3d x_B = Vector3d(0.5, 0.8, 0) + 0.001 * zeta;
	Vector3d eta_B = this->SP.to_eta2(x_B);
	Vector3d eta_B_ex = Vector3d(Unit::deg2rad(300), 0, 0.001);
	EXPECT_NEAR(eta_B.x(), eta_B_ex.x(), 1e-3);
	EXPECT_NEAR(eta_B.y(), eta_B_ex.y(), 1e-7);
	EXPECT_NEAR(eta_B.z(), eta_B_ex.z(), 1e-7);

	// 点B2(ζ方向に10mm移動，精度が少し落ちる)
	Vector3d x_B2 = Vector3d(0.5, 0.8, 0) + 0.01 * zeta;
	Vector3d eta_B2 = this->SP.to_eta2(x_B2);
	Vector3d eta_B2_ex = Vector3d(Unit::deg2rad(300), 0, 0.01);
	EXPECT_NEAR(eta_B2.x(), eta_B2_ex.x(), 1e-3);
	EXPECT_NEAR(eta_B2.y(), eta_B2_ex.y(), 1e-6);
	EXPECT_NEAR(eta_B2.z(), eta_B2_ex.z(), 1e-6);

	// 点D(η方向にr/2だけ移動)
	Vector3d x_D = Vector3d(0.5, 0.4, 0);
	Vector3d eta_D = this->SP.to_eta2(x_D);
	Vector3d eta_D_ex = Vector3d(Unit::deg2rad(300), 0.4, 0);
	EXPECT_NEAR(eta_D.x(), eta_D_ex.x(), 1e-3);
	EXPECT_NEAR(eta_D.y(), eta_D_ex.y(), 1e-7);
	EXPECT_NEAR(eta_D.z(), eta_D_ex.z(), 1e-7);
	return;
}

// (1)-2 "座標"を直交座標系⇒螺旋座標系に変換，手計算と比較
TEST_F(SpiralTest, x_to_eta_02) {
	this->init_data02();
	// 点A(x=0から60°分移動した位置)
	Vector3d x_A = Vector3d(0.01, 0.05, 0.0866025403);
	Vector3d eta_A = this->SP.to_eta2(x_A);
	Vector3d eta_A_ex = Vector3d(Unit::deg2rad(60), 0, 0);
	EXPECT_NEAR(eta_A.x(), eta_A_ex.x(), 1e-3);
	EXPECT_NEAR(eta_A.y(), eta_A_ex.y(), 1e-7);
	EXPECT_NEAR(eta_A.z(), eta_A_ex.z(), 1e-7);

	Vector3d zeta = Vector3d(0.995471495, 0.08232483, -0.047530263); // ζの単位ベクトル

	// 点B(ζ方向に0.1mm移動)
	Vector3d x_B = Vector3d(0.01, 0.05, 0.0866025403) + 0.0001 * zeta;
	Vector3d eta_B = this->SP.to_eta2(x_B);
	Vector3d eta_B_ex = Vector3d(Unit::deg2rad(60), 0, 0.0001);
	EXPECT_NEAR(eta_B.x(), eta_B_ex.x(), 1e-3);
	EXPECT_NEAR(eta_B.y(), eta_B_ex.y(), 1e-7);
	EXPECT_NEAR(eta_B.z(), eta_B_ex.z(), 1e-7);

	// 点B2(ζ方向に1mm移動，精度が少し落ちる)
	Vector3d x_B2 = Vector3d(0.01, 0.05, 0.0866025403) + 0.001 * zeta;
	Vector3d eta_B2 = this->SP.to_eta2(x_B2);
	Vector3d eta_B2_ex = Vector3d(Unit::deg2rad(60), 0, 0.001);
	EXPECT_NEAR(eta_B2.x(), eta_B2_ex.x(), 1e-3);
	EXPECT_NEAR(eta_B2.y(), eta_B2_ex.y(), 1e-6);
	EXPECT_NEAR(eta_B2.z(), eta_B2_ex.z(), 1e-6);

	// 点D(η方向にr/2だけ移動)
	Vector3d x_D = Vector3d(0.01, 0.05 / 2, 0.0866025403 / 2);
	Vector3d eta_D = this->SP.to_eta2(x_D);
	Vector3d eta_D_ex = Vector3d(Unit::deg2rad(60), 0.05, 0);
	EXPECT_NEAR(eta_D.x(), eta_D_ex.x(), 1e-3);
	EXPECT_NEAR(eta_D.y(), eta_D_ex.y(), 1e-7);
	EXPECT_NEAR(eta_D.z(), eta_D_ex.z(), 1e-7);
	return;
}

// (2) 螺旋座標に変換した"座標"をもとの直交座標系に変換できるか
TEST_F(SpiralTest, x_to_eta_to_x) {
	this->init_data01();
	Vector3d zeta = Vector3d(0.992951092, 0, -0.118524806); // ζの単位ベクトル

	// 点B(ζ方向に1mm移動)
	Vector3d x_B = Vector3d(0.5, 0.8, 0) + 0.001 * zeta;
	Vector3d eta_B = this->SP.to_eta2(x_B);
	Vector3d x_B_ = this->SP.to_xyz(eta_B);
	EXPECT_LT((x_B - x_B_).norm(), 1e-7);

	// 点D(η方向にr/2だけ移動)
	Vector3d x_D = Vector3d(0.5, 0.4, 0);
	Vector3d eta_D = this->SP.to_eta2(x_D);
	Vector3d x_D_ = this->SP.to_xyz(eta_D);
	EXPECT_LT((x_D - x_D_).norm(), 1e-7);

	return;
};

// (3)-1 "ベクトル"を直交座標系⇒螺旋座標系に変換
TEST_F(SpiralTest, get_xyz2eta_01) {
	this->init_data01();

	// (1)-1の点Aを基準にした螺旋座標系で考える
	// eta, zeta, xaiの単位ベクトル（直交座標系）を変換できるか確認
	double theta = Unit::deg2rad(300);
	Vector3d eta  = Vector3d(0, -1, 0);	// ηの単位ベクトル
	Vector3d zeta = Vector3d(0.992951092, 0, -0.118524806); // ζの単位ベクトル
	Vector3d xi   = eta.cross(zeta);	// ξの単位ベクトル

	Matrix3d xyz2eta = this->SP.get_xyz2eta(theta);
	Vector3d eta_ = xyz2eta * eta;
	Vector3d zeta_ = xyz2eta * zeta;
	Vector3d xi_ = xyz2eta * xi;
	Vector3d xi_ex(1, 0, 0), eta_ex(0, 1, 0), zeta_ex(0, 0, 1);
	EXPECT_NEAR((xi_ - xi_ex).norm(), 0, 1e-7);
	EXPECT_NEAR((eta_ - eta_ex).norm(), 0, 1e-7);
	EXPECT_NEAR((zeta_ - zeta_ex).norm(), 0, 1e-7);

	return;
}

// (3)-2 "ベクトル"を直交座標系⇒螺旋座標系に変換
TEST_F(SpiralTest, get_xyz2eta_02) {
	this->init_data02();

	// (1)-1の点Aを基準にした螺旋座標系で考える
	// eta, zeta, xaiの単位ベクトル（直交座標系）を変換できるか確認
	double theta = Unit::deg2rad(60);
	Vector3d eta = Vector3d(0, -0.5, -0.866025403);	// ηの単位ベクトル
	Vector3d zeta = Vector3d(0.995471495, 0.08232483, -0.047530263); // ζの単位ベクトル
	Vector3d xi = eta.cross(zeta);	// ξの単位ベクトル

	Matrix3d xyz2eta = this->SP.get_xyz2eta(theta);
	Vector3d eta_ = xyz2eta * eta;
	Vector3d zeta_ = xyz2eta * zeta;
	Vector3d xi_ = xyz2eta * xi;
	Vector3d xi_ex(1, 0, 0), eta_ex(0, 1, 0), zeta_ex(0, 0, 1);
	EXPECT_NEAR((xi_ - xi_ex).norm(), 0, 1e-7);
	EXPECT_NEAR((eta_ - eta_ex).norm(), 0, 1e-7);
	EXPECT_NEAR((zeta_ - zeta_ex).norm(), 0, 1e-7);


	Matrix3d eta2xyz = xyz2eta.inverse();
	Vector3d zeta__ = eta2xyz * zeta_ex;
	Vector3d xi__ = eta2xyz * xi_ex;
	return;
}



// (4) (3)の変換行列が回転行列の性質を満たしているか判定
TEST_F(SpiralTest, get_xyz2eta_10) {
	this->init_data01();
	Vector3d x = Vector3d::Ones();

	Vector3d eta = this->SP.to_eta2(x);

	Matrix3d xyz2eta = this->SP.get_xyz2eta(eta[0]);

	EXPECT_GT(1e-2, xyz2eta.determinant() - 1.0);		// 行列式が１となることの確認．

	Matrix3d inverse = xyz2eta.inverse();
	Matrix3d transpose = xyz2eta.transpose();

	EXPECT_GT(1e-2, (inverse - transpose).norm());		// 逆行列と転地の一致の確認．

	return;
}

// getterのテスト
TEST_F(SpiralTest, getter) {
	this->init_data01();
	double nd = this->SP.get_nd();
	EXPECT_GT(1e-2, nd - 1.0);		// 手計算との一致の確認．

	double r = this->SP.get_r();
	EXPECT_GT(1e-2, r - 0.8);		// 設定値との一致の確認．

	return;
};


