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

// ���W�𗆐����W�n�ɕϊ����Ė߂��i�ڍׂȃe�X�g��Spiral�N���X�Ŏ��{�ς݂ł��邽�߁C�P���Ȃ��̂̂ݎ��{�j
TEST_F(BS_CylinderTest, to_etacoord0) {
	this->init();

	double th = 1e1;

	Vector3d eta0(th, 1e-3, 2e-3);
	Vector3d x0 = this->ST.CY->to_inertialcoord(0, eta0);

	Vector3d eta1 = this->ST.CY->to_etacoord(0, x0);

	double delta = (eta0 - eta1).norm();

	EXPECT_GT(1e-2, delta);
};

// ���x�x�N�g���𗆐����W�n�ɕϊ����Ė߂�
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

// �ʑ��p0���̈ʒu�ŋʂƍa���ڐG�����ۂ̐ڐG�ʒu�𓱏o
TEST_F(BS_CylinderTest, calc_slice_test1) {
	this->init();

	int is = 0, ig = 0; // ��ԍ��C�a�ԍ�
	double th = 0;	//	�����ʑ��p
	Vector3d p(0, 0.8 - 1e-6, 0);	// �ڐG�_[m]
	Vector3d t(0.8, 0, 2 * Numeric::pi * 1.0);	// �����i�s�x�N�g��[m]
	Vector3d xai = t / t.norm();				// �����i�s�����P�ʃx�N�g��[-]
	double a = 1e-4;							// �ڐG�ȉ~���a�i�K���j
	int ms = 21;								// �X���C�X������							
	double bl_r = 0.15;							// �ʔ��a[m](�K��)
	Vector3d ps[21];

	this->ST.CY->calc_slice(is, ig, th, p, xai, a, ms, bl_r, ps);

	// (1) �ڐG�ȉ~�̒[����[�ƁC�g�[���X���S����ڐG�_�ւ̃x�N�g�������s���Ă��邩�̂��̊m�F
	Vector3d op = Vector3d(0, 1, 0);
	double ip1 = op.normalized().dot((ps[20] - ps[0]).normalized());
	EXPECT_NEAR(ip1, 0, 1e-6);

	// (2) �ڐG�ȉ~�̒[����[�ƁC�����i�s���������s���Ă��邩�̊m�F
	double ip2 = xai.dot((ps[20] - ps[0]).normalized());
	EXPECT_NEAR(ip2, 0, 1e-6);

	// (3) �ڐG�ȉ~�̒[����[�̋������C�ڐG�ȉ~�̒����ɓ��������m�F
	// (�X���C�X��������Z���Ȃ邽�߁C�������Čv�Z)
	double ip3 = (ps[20] - ps[0]).norm();
	EXPECT_NEAR(ip3, a * 2 / 21.0 * 20.0 , 1e-6);

	// (4) �ڐG�ȉ~�̕В[����aR���S�ւ̋������������������̊m�F
	Vector3d Og(0, 1.0, 0);	// �aR���S[m]
	double l1 = (ps[20] - Og).norm();
	double l2 = (ps[0] - Og).norm();
	EXPECT_NEAR(l1, l2, 1e-6);

	// (5) �ڐG�ȉ~�̕В[����aR���S�ւ̋������C�ڐG�_���璆�S�ւ̋����̕����傫�����̊m�F
	double l0 = (ps[20] - Og).norm();
	EXPECT_GE(l0, l1);

	// (6) �ڐG�ȉ~�̕В[����aR���S�ւ̋������C�U�b�N���a���a�Ɠ��������̊m�F
	EXPECT_NEAR(l1, 0.2, 1e-5);

	return;
}
