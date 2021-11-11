/*******************************************************************************
!						"BS_BallScrewTest.cpp"
!														2021/05/19	[Core-T] ����
!	�����l�v�Z�̃��\�b�h�����C���Ƀe�X�g�D���́E�o�͗p�̊֐��͐������@�\������̂Ƃ��čl����D
!	�������C�ʏ������W�Ȃǎ኱���G�Ȍv�Z�����Ă�����̂̓e�X�g����D
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

	// �{�[���˂��̊􉽌`����`
	// �a��C�����̌X����30���̒P���`��ɂ���
	void SingleNut_init() {
		// �������_
		IN.tribology.rollingresistance = BS_In::Tribology::RollingResistance::RollingResistanceNothing;
		IN.tribology.coulomb = BS_In::Tribology::Tangent;
		IN.tribology.filmThickness = BS_In::Tribology::FilmThicknessNothing;
		IN.tribology.coulomb_slope = 1000;

		// ��H�E��
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

		// �V���t�g�����l(�����̏����ʑ��p��0��)
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
		// �aR���SPCD[m]�i�apcd���a�C0.0357535 m ���炢�j
		this->IN.shaft[0].spiral[0].groove[0].eta[0] = -0.5e-3;
		this->IN.shaft[0].spiral[0].groove[1].eta[0] = -0.5e-3;

		this->IN.shaft[0].spiral[0].groove[0].sigma = 0.00000006;
		this->IN.shaft[0].spiral[0].groove[1].sigma = 0.00000006;
		this->IN.shaft[0].spiral[0].alp = 0;
		this->IN.shaft[0].spiral[0].r = 50e-3;
		this->IN.shaft[0].spiral[0].l = 181.2879e-3;// (= 2��r * tan30��)

		// �i�b�g�����l
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
		this->IN.nut[0].spiral[0].l = 181.2879e-3;// (= 2��r * tan30��)

		this->IN.nut[0].m = 1;
		this->IN.nut[0].young = 207760000000;
		this->IN.nut[0].poisson = 0.3;
		this->IN.nut[0].Ix = 1, this->IN.nut[0].Iyz = 1, this->IN.nut[0].Iyz = 1;

		// Pair�N���X�����l
		int n_ball = this->IN.circuit[0].ball.size();
		this->IN.BallShaftPair.resize(n_ball);
		for (int i = 0; i < n_ball; i++) {
			this->IN.BallShaftPair[0].groove[0].zeta = 0.2;	// ������
			this->IN.BallShaftPair[0].groove[0].mu = 0.1;		// ���C�W��
			this->IN.BallShaftPair[0].groove[1].zeta = 0.2;	// ������
			this->IN.BallShaftPair[0].groove[1].mu = 0.1;		// ���C�W��
		}
		this->IN.BallNutPair.resize(n_ball);
		for (int i = 0; i < n_ball; i++) {
			this->IN.BallNutPair[0].groove[0].zeta = 0.2;	// ������
			this->IN.BallNutPair[0].groove[0].mu = 0.1;		// ���C�W��
			this->IN.BallNutPair[0].groove[1].zeta = 0.2;	// ������
			this->IN.BallNutPair[0].groove[1].mu = 0.1;		// ���C�W��
		}

		this->IN.oil.eta = 0.5;	// �S�x�iPa*s�j(�K��)
		this->IN.oil.beta = 0.5;	// ���x�S�x�W���i�K���j
		this->IN.oil.k = 0.145;	// ���M�`����[W/(m�EK)]
		this->IN.oil.alpha = 0.2;	// ���͔S�x�W��[mm2/kgf]
		this->IN.oil.lm = -1;	// ���j�X�J�X����

		// �׏d�i�K�v�ɉ����Ċe�e�X�g�P�[�X�̒��Œ�`�j
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

		// �������_
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

// (�ʏ����l�v�Z) �ʂ����Ԋu�ɕ���ł��邩�m�F
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

// (�����l�v�Z�P�[�X0) �׏d0�̂Ƃ��͉������Ȃ��D
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
// (�����l�v�Z�P�[�X1) ���A�L�V�A���̂Ƃ��C��̃A�L�V�A�������ܕ��V���t�g�Ƌʂ��ړ����邩�m�F
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

	double l = 181.2879e-3;		// ���[�h��[m]
	for (int i = 0; i < 5; i++) {
		double x = this->OUT.CC[0].BL[i].x[0];
		double y = this->OUT.CC[0].BL[i].x[1];
		double z = this->OUT.CC[0].BL[i].x[2];
		double _x = l / 4 * i + 0.86602540e-3;
		double _th = Numeric::pi * 2 / 4 * i;
		double _y = cos(_th) * 50.5e-3;
		double _z = sin(_th) * 50.5e-3;
		// ���e�덷��y, z���W��1%���x
		EXPECT_NEAR(this->OUT.CC[0].BL[i].x[0], _x, 5e-4);
		EXPECT_NEAR(this->OUT.CC[0].BL[i].x[1], _y, 5e-4);
		EXPECT_NEAR(this->OUT.CC[0].BL[i].x[2], _z, 5e-4);
	}
	return;
}

// (�����l�v�Z�P�[�X2) �����W�A���̂Ƃ��C��̃��W�A�������ܕ��V���t�g���ړ����邩�m�F
// �܂��C���ׂĂ̋ʂ��i�b�g�ɐڐG���Ă���C���Ԋu�ɂȂ��Ă��邩�m�F
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

	double l = 181.2879e-3;		// ���[�h��[m]
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

// (�����l�v�Z�P�[�X3) �~�X�A���C�����g������Ƃ��C�׏d�������������Ɏ����X�����m�F
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