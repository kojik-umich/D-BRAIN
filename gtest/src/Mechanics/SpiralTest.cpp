#include "pch.h"
/*******************************************************************************
!								"SpiralTest.cpp"
!												2021/4/2	[Core-T]	���C����
!	�������W�n�̗������W�n�̕ϊ����e�X�g���Ă���D
!	�e�X�g���ڂ�
!	(1) "���W"�𒼌����W�n�˗������W�n�ɕϊ�
!	(2) �������W�ɕϊ�����"���W"�����Ƃ̒������W�n�ɕϊ��ł��邩
!	(3) "�x�N�g��"�𒼌����W�n�˗������W�n�ɕϊ�
!	(4) (3)�̕ϊ��s�񂪉�]�s��̐����𖞂����Ă��邩����
!	
!	��B:\�\�[�X�R�[�h��\D-BRAIN\docs\26_D-BS\03_�]��\01_�P�̃e�X�g\SpiralTest.pptx�ɕ⑫�������L��
!*******************************************************************************/


class SpiralTest : public ::testing::Test {
public:
	Spiral SP;
	
	virtual void SetUp() {
	}

	// �e�X�g�f�[�^1
	void init_data01() {
		double alp = Unit::deg2rad(60);// �����ʑ��p[rad](30��)
		double l = 0.6;
		double r = 0.8;

		this->SP.init(alp, l, r);
	}

	// �e�X�g�f�[�^2
	void init_data02() {
		double alp = Unit::deg2rad(0);	// �����ʑ��p[rad](0��)
		double l = 0.06;
		double r = 0.1;

		this->SP.init(alp, l, r);
	}
};
// (1)-1 "���W"�𒼌����W�n�˗������W�n�ɕϊ��C��v�Z�Ɣ�r
TEST_F(SpiralTest, x_to_eta_01) {
	this->init_data01();
	// �_A(x=0����300�����ړ������ʒu)
	Vector3d x_A = Vector3d(0.5, 0.8 , 0);
	Vector3d eta_A = this->SP.to_eta2(x_A);
	Vector3d eta_A_ex = Vector3d(Unit::deg2rad(300), 0, 0);
	EXPECT_NEAR(eta_A.x(), eta_A_ex.x(), 1e-3);
	EXPECT_NEAR(eta_A.y(), eta_A_ex.y(), 1e-7);
	EXPECT_NEAR(eta_A.z(), eta_A_ex.z(), 1e-7);
	Vector3d x1 = this->SP.to_xyz(eta_A);
	Vector3d x2 = this->SP.to_xyz(eta_A_ex);
	Vector3d zeta = Vector3d(0.992951092, 0, -0.118524806); // �Ă̒P�ʃx�N�g��

	// �_B(�ĕ�����1mm�ړ�)
	Vector3d x_B = Vector3d(0.5, 0.8, 0) + 0.001 * zeta;
	Vector3d eta_B = this->SP.to_eta2(x_B);
	Vector3d eta_B_ex = Vector3d(Unit::deg2rad(300), 0, 0.001);
	EXPECT_NEAR(eta_B.x(), eta_B_ex.x(), 1e-3);
	EXPECT_NEAR(eta_B.y(), eta_B_ex.y(), 1e-7);
	EXPECT_NEAR(eta_B.z(), eta_B_ex.z(), 1e-7);

	// �_B2(�ĕ�����10mm�ړ��C���x������������)
	Vector3d x_B2 = Vector3d(0.5, 0.8, 0) + 0.01 * zeta;
	Vector3d eta_B2 = this->SP.to_eta2(x_B2);
	Vector3d eta_B2_ex = Vector3d(Unit::deg2rad(300), 0, 0.01);
	EXPECT_NEAR(eta_B2.x(), eta_B2_ex.x(), 1e-3);
	EXPECT_NEAR(eta_B2.y(), eta_B2_ex.y(), 1e-6);
	EXPECT_NEAR(eta_B2.z(), eta_B2_ex.z(), 1e-6);

	// �_D(�ŕ�����r/2�����ړ�)
	Vector3d x_D = Vector3d(0.5, 0.4, 0);
	Vector3d eta_D = this->SP.to_eta2(x_D);
	Vector3d eta_D_ex = Vector3d(Unit::deg2rad(300), 0.4, 0);
	EXPECT_NEAR(eta_D.x(), eta_D_ex.x(), 1e-3);
	EXPECT_NEAR(eta_D.y(), eta_D_ex.y(), 1e-7);
	EXPECT_NEAR(eta_D.z(), eta_D_ex.z(), 1e-7);
	return;
}

// (1)-2 "���W"�𒼌����W�n�˗������W�n�ɕϊ��C��v�Z�Ɣ�r
TEST_F(SpiralTest, x_to_eta_02) {
	this->init_data02();
	// �_A(x=0����60�����ړ������ʒu)
	Vector3d x_A = Vector3d(0.01, 0.05, 0.0866025403);
	Vector3d eta_A = this->SP.to_eta2(x_A);
	Vector3d eta_A_ex = Vector3d(Unit::deg2rad(60), 0, 0);
	EXPECT_NEAR(eta_A.x(), eta_A_ex.x(), 1e-3);
	EXPECT_NEAR(eta_A.y(), eta_A_ex.y(), 1e-7);
	EXPECT_NEAR(eta_A.z(), eta_A_ex.z(), 1e-7);

	Vector3d zeta = Vector3d(0.995471495, 0.08232483, -0.047530263); // �Ă̒P�ʃx�N�g��

	// �_B(�ĕ�����0.1mm�ړ�)
	Vector3d x_B = Vector3d(0.01, 0.05, 0.0866025403) + 0.0001 * zeta;
	Vector3d eta_B = this->SP.to_eta2(x_B);
	Vector3d eta_B_ex = Vector3d(Unit::deg2rad(60), 0, 0.0001);
	EXPECT_NEAR(eta_B.x(), eta_B_ex.x(), 1e-3);
	EXPECT_NEAR(eta_B.y(), eta_B_ex.y(), 1e-7);
	EXPECT_NEAR(eta_B.z(), eta_B_ex.z(), 1e-7);

	// �_B2(�ĕ�����1mm�ړ��C���x������������)
	Vector3d x_B2 = Vector3d(0.01, 0.05, 0.0866025403) + 0.001 * zeta;
	Vector3d eta_B2 = this->SP.to_eta2(x_B2);
	Vector3d eta_B2_ex = Vector3d(Unit::deg2rad(60), 0, 0.001);
	EXPECT_NEAR(eta_B2.x(), eta_B2_ex.x(), 1e-3);
	EXPECT_NEAR(eta_B2.y(), eta_B2_ex.y(), 1e-6);
	EXPECT_NEAR(eta_B2.z(), eta_B2_ex.z(), 1e-6);

	// �_D(�ŕ�����r/2�����ړ�)
	Vector3d x_D = Vector3d(0.01, 0.05 / 2, 0.0866025403 / 2);
	Vector3d eta_D = this->SP.to_eta2(x_D);
	Vector3d eta_D_ex = Vector3d(Unit::deg2rad(60), 0.05, 0);
	EXPECT_NEAR(eta_D.x(), eta_D_ex.x(), 1e-3);
	EXPECT_NEAR(eta_D.y(), eta_D_ex.y(), 1e-7);
	EXPECT_NEAR(eta_D.z(), eta_D_ex.z(), 1e-7);
	return;
}

// (2) �������W�ɕϊ�����"���W"�����Ƃ̒������W�n�ɕϊ��ł��邩
TEST_F(SpiralTest, x_to_eta_to_x) {
	this->init_data01();
	Vector3d zeta = Vector3d(0.992951092, 0, -0.118524806); // �Ă̒P�ʃx�N�g��

	// �_B(�ĕ�����1mm�ړ�)
	Vector3d x_B = Vector3d(0.5, 0.8, 0) + 0.001 * zeta;
	Vector3d eta_B = this->SP.to_eta2(x_B);
	Vector3d x_B_ = this->SP.to_xyz(eta_B);
	EXPECT_LT((x_B - x_B_).norm(), 1e-7);

	// �_D(�ŕ�����r/2�����ړ�)
	Vector3d x_D = Vector3d(0.5, 0.4, 0);
	Vector3d eta_D = this->SP.to_eta2(x_D);
	Vector3d x_D_ = this->SP.to_xyz(eta_D);
	EXPECT_LT((x_D - x_D_).norm(), 1e-7);

	return;
};

// (3)-1 "�x�N�g��"�𒼌����W�n�˗������W�n�ɕϊ�
TEST_F(SpiralTest, get_xyz2eta_01) {
	this->init_data01();

	// (1)-1�̓_A����ɂ����������W�n�ōl����
	// eta, zeta, xai�̒P�ʃx�N�g���i�������W�n�j��ϊ��ł��邩�m�F
	double theta = Unit::deg2rad(300);
	Vector3d eta  = Vector3d(0, -1, 0);	// �ł̒P�ʃx�N�g��
	Vector3d zeta = Vector3d(0.992951092, 0, -0.118524806); // �Ă̒P�ʃx�N�g��
	Vector3d xi   = eta.cross(zeta);	// �̂̒P�ʃx�N�g��

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

// (3)-2 "�x�N�g��"�𒼌����W�n�˗������W�n�ɕϊ�
TEST_F(SpiralTest, get_xyz2eta_02) {
	this->init_data02();

	// (1)-1�̓_A����ɂ����������W�n�ōl����
	// eta, zeta, xai�̒P�ʃx�N�g���i�������W�n�j��ϊ��ł��邩�m�F
	double theta = Unit::deg2rad(60);
	Vector3d eta = Vector3d(0, -0.5, -0.866025403);	// �ł̒P�ʃx�N�g��
	Vector3d zeta = Vector3d(0.995471495, 0.08232483, -0.047530263); // �Ă̒P�ʃx�N�g��
	Vector3d xi = eta.cross(zeta);	// �̂̒P�ʃx�N�g��

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



// (4) (3)�̕ϊ��s�񂪉�]�s��̐����𖞂����Ă��邩����
TEST_F(SpiralTest, get_xyz2eta_10) {
	this->init_data01();
	Vector3d x = Vector3d::Ones();

	Vector3d eta = this->SP.to_eta2(x);

	Matrix3d xyz2eta = this->SP.get_xyz2eta(eta[0]);

	EXPECT_GT(1e-2, xyz2eta.determinant() - 1.0);		// �s�񎮂��P�ƂȂ邱�Ƃ̊m�F�D

	Matrix3d inverse = xyz2eta.inverse();
	Matrix3d transpose = xyz2eta.transpose();

	EXPECT_GT(1e-2, (inverse - transpose).norm());		// �t�s��Ɠ]�n�̈�v�̊m�F�D

	return;
}

// getter�̃e�X�g
TEST_F(SpiralTest, getter) {
	this->init_data01();
	double nd = this->SP.get_nd();
	EXPECT_GT(1e-2, nd - 1.0);		// ��v�Z�Ƃ̈�v�̊m�F�D

	double r = this->SP.get_r();
	EXPECT_GT(1e-2, r - 0.8);		// �ݒ�l�Ƃ̈�v�̊m�F�D

	return;
};


