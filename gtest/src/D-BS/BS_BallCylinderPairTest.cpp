/*******************************************************************************
!						"BS_BallCylinderPairTest.cpp"
!														2021/05/13	[Core-T] ����
!	��{�I�� D-B4P �� B4P_BallRingPairTest.cpp �̃e�X�g�P�[�X�𗬗p
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

	// �{�[���˂��̊􉽌`����`
	// �a��C�����̌X����30���̒P���`��ɂ���
	void BS_init() {
		// Pair�N���X�����l
		this->IN.BallShaftPair.resize(1);
		this->IN.BallShaftPair[0].groove[0].zeta = 0.2;	// ������
		this->IN.BallShaftPair[0].groove[0].mu = 0.1;		// ���C�W��
		this->IN.BallShaftPair[0].groove[1].zeta = 0.2;	// ������
		this->IN.BallShaftPair[0].groove[1].mu = 0.1;		// ���C�W��
		// �ʕ����l
		double D = 20e-3;				// �ʌa[m]
		double E = 207760000000;		// �����O��[Pa]
		double por = 0.3;				// �|�A�\����[-]
		double den = 7830;				// ���x[kg/m^3]
		double rms = 0.00000002;		// �e��rms[m]
		this->BL = new Ball;
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}
		this->BL->init(D, E, por, den, rms, x_const, Rx_const);

		// �V���t�g�����l
		IN.shaft.resize(1);
		IN.shaft[0].spiral.resize(1);
		IN.shaft[0].spiral[0].groove[0].eta[1] = 0;
		IN.shaft[0].spiral[0].groove[1].eta[1] = 0;
		IN.shaft[0].spiral[0].groove[0].r = 11e-3;
		IN.shaft[0].spiral[0].groove[1].r = 1;
		// �aR���SPCD[m]�i�apcd���a�C0.0357535 m ���炢�j
		IN.shaft[0].spiral[0].groove[0].eta[0] = 0;
		IN.shaft[0].spiral[0].groove[1].eta[0] = 0;

		IN.shaft[0].spiral[0].groove[0].sigma = 0.00000006;
		IN.shaft[0].spiral[0].groove[1].sigma = 0.00000006;
		IN.shaft[0].spiral[0].alp = 0.5 * Numeric::pi;
		IN.shaft[0].spiral[0].r = 50e-3;
		IN.shaft[0].spiral[0].l = 181.2879e-3;// (= 2��r * tan30��)

		IN.shaft[0].m = 1;
		IN.shaft[0].young = 207760000000;
		IN.shaft[0].poisson = 0.3;
		IN.shaft[0].Ix = 1, IN.shaft[0].Iyz = 1, IN.shaft[0].Iyz = 1;

		this->IN.oil.eta = 0.5;	// �S�x�iPa*s�j(�K��)
		this->IN.oil.beta = 0.5;	// ���x�S�x�W���i�K���j
		this->IN.oil.k = 0.145;	// ���M�`����[W/(m�EK)]
		this->IN.oil.alpha = 0.2;	// ���͔S�x�W��[mm2/kgf]
		this->IN.oil.lm = -1;	// ���j�X�J�X����


		// �������_
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


		// �O�֕����l
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
		IN.nut[0].spiral[0].l = 181.2879e-3;// (= 2��r * tan30��)

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
		//this->BNP.TB.FT = new Tribology::FilmThicknessNothing();	// ����������0�ɌŒ�
		//this->BNP.TB.cs = 1000;
	}

	// 4�ʂ̎��󏏌��i25BSWZ01�j�𗬗p
	void init_25BSWZ() {
		// ���̓f�[�^�͉��L�̂��̂�p����
		// \\ans00978\kiken\BRAIN\�\�[�X�R�[�h��\DBRAINpp\docs\06_D4Bv100\05_�P�̃e�X�g\B4P_BallRingPair
		// ����^�ԁ@25BSWZ01
		// Pair�N���X�����l
		this->IN.BallShaftPair.resize(1);
		this->IN.BallShaftPair[0].groove[0].zeta = 0.2;	// ������
		this->IN.BallShaftPair[0].groove[0].mu = 0.1;		// ���C�W��
		this->IN.BallShaftPair[0].groove[1].zeta = 0.2;	// ������
		this->IN.BallShaftPair[0].groove[1].mu = 0.1;		// ���C�W��
		// �ʕ����l
		double D = 0.00635;				// �ʌa[m]
		double E = 207760000000;		// �����O��[Pa]
		double por = 0.29;				// �|�A�\����[-]
		double den = 7830;				// ���x[kg/m^3]
		double rms = 0.00000002;		// �e��rms[m]
		this->BL = new Ball;
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}
		this->BL->init(D, E, por, den, rms, x_const, Rx_const);

		// ���֕����l
		IN.shaft.resize(1);
		IN.shaft[0].spiral.resize(1);
		IN.shaft[0].spiral[0].r = 0.0355 / 2;
		IN.shaft[0].spiral[0].groove[0].eta[1] = 0.0000725;
		IN.shaft[0].spiral[0].groove[1].eta[1] = -0.0000725;
		IN.shaft[0].spiral[0].groove[0].r = 0.00332105;
		IN.shaft[0].spiral[0].groove[1].r = 0.00332105;
		// �aR���SPCD[m]�i�apcd���a�C0.0357535 m ���炢�j
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

		this->IN.oil.eta = 0.5;	// �S�x�iPa*s�j(�K��)
		this->IN.oil.beta = 0.5;	// ���x�S�x�W���i�K���j
		this->IN.oil.k = 0.145;	// ���M�`����[W/(m�EK)]
		this->IN.oil.alpha = 0.2;	// ���͔S�x�W��[mm2/kgf]
		this->IN.oil.lm = -1;	// ���j�X�J�X����


		// �������_
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


		// �O�֕����l
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

// (�׏d�v�Z�P�[�X1) �V���v���ȃP�[�X��z��
// �ڐG�p180�����ʈʑ��p0���Ŗ���]�ŕǖʂɑ΂��ĕ��i�^�����Ă��鎞������
TEST_F(BS_BallCylinderPairTest, calc_force_case1) {
	this->BS_init();
	// �ʐڋߗ� 1um�C�ʂ̐i�s�����͗����ڐ������ɂȂ�悤�ɐݒ�
	Vector3d bl_x = Vector3d(0, 0, 50e-3 + 1.001e-3);
	Vector3d bl_v = Vector3d(0.5, 0.8660254, 0);
	Vector3d bl_w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // �����N�H�[�^�j�I���͑S�� [0,0,0,1] �Ƃ���D�i�������R���X�g���N�^�̎d�l�㏇�Ԃ��قȂ�j

	this->BL->set_y(bl_x, bl_v, q, bl_w);


	Vector3d Fbi, Tbi, Fib, Tib, Fs, Ts;

	this->BNP.get_FT(Fbi, Tbi, Tib, Fs, Ts);

	// 1. �ʂ���i�b�g�ɓ����w���c�׏d�� -�� �����i= +z �����j�̐��������ɂȂ��Ă��邩�m�F
	Vector3d Fn = -Fbi - Fs;
	double err = 1e-3;
	EXPECT_NEAR(Fn.z(), Fn.norm(), Fn.norm() * err);
	// 2. ���C�� = �w���c�׏d �~ ���C�W�� �ɂȂ��Ă��邩�m�F
	double F = Fs.norm();
	EXPECT_NEAR(Fs.norm(), Fn.norm() * 0.1, Fs.norm() * err);
	// 3. �ʂ���i�b�g�ɓ������C�͋ʂ��i�s��������Ƃ͓��������ɂȂ��Ă��邱�Ƃ��m�F
	EXPECT_NEAR((Fs.normalized() - bl_v.normalized()).norm(), 0, err);


	return;
}

// (�׏d�v�Z�P�[�X2) �ʑ��p�E�ڐG�p����[���̃P�[�X������
// �ڐG�p150�����ʈʑ��p90���Ŗ���]�ŕǖʂɑ΂��ĕ��i�^�����Ă��鎞������
TEST_F(BS_BallCylinderPairTest, calc_force_case2) {
	this->BS_init();
	// �ʐڋߗ� 1um�C�ʂ̐i�s�����͗����ڐ������ɂȂ�悤�ɐݒ�
	// ���̗����̏����ʑ��p��90���Ȃ̂ŁC�����猩��-y�����ɋʂ�����
	Vector3d bl_x = Vector3d(181.2879e-3 * 0.25 + 1.001e-3 * sin30 * cos30, -50e-3 - 1.001e-3 * cos30, 1.001e-3 * sin30 * sin30);
	Vector3d bl_v = Vector3d(sin30, 0, -cos30);
	Vector3d bl_w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // �����N�H�[�^�j�I���͑S�� [0,0,0,1] �Ƃ���D�i�������R���X�g���N�^�̎d�l�㏇�Ԃ��قȂ�j

	this->BL->set_y(bl_x, bl_v, q, bl_w);


	Vector3d Fbi, Tbi, Fib, Tib, Fs, Ts;

	this->BNP.get_FT(Fbi, Tbi, Tib, Fs, Ts);

	// 1. �ʂ���i�b�g�ɓ����w���c�׏d�̓łɑ΂���150���̌����ɂȂ��Ă��邩�m�F
	Vector3d e = (-Fbi - Fs).normalized();
	double err = 1e-3;
	Vector3d e_ex = Vector3d(sin30 * cos30, -cos30, sin30 * sin30);
	IOFormat CleanFmt(4, 0, ", ", ",", "", "", "[", "]");
	EXPECT_NEAR((e_ex - e).norm(), 0, err) << "e = " << e.format(CleanFmt) << endl << "e_ex = " << e_ex.format(CleanFmt);
	// 2. ���C�� = �w���c�׏d �~ ���C�W�� �ɂȂ��Ă��邩�m�F
	Vector3d Fn = -Fbi - Fs;
	EXPECT_NEAR(Fs.norm(), Fn.norm() * 0.1, Fs.norm() * err);
	// 3. �ʂ���i�b�g�ɓ������C�͋ʂ��i�s��������Ƃ͓��������ɂȂ��Ă��邱�Ƃ��m�F
	EXPECT_NEAR((Fs.normalized() - bl_v.normalized()).norm(), 0, err);

	return;
}

// (�׏d�v�Z�P�[�X3) �ʂ��X�s�������ɉ�]���Ă���P�[�X������
// �ڐG�p0�����ʈʑ��p0���ŋʂ��ǖʂɑ΂��ăX�s�����Ă��鎞������
TEST_F(BS_BallCylinderPairTest, calc_force_case3) {
	this->BS_init();
	// �ʐڋߗ� 1um�C�ʂ̐i�s�����͗����ڐ������ɂȂ�悤�ɐݒ�
	// ���̗����̏����ʑ��p��90���Ȃ̂ŁC�����猩��+z�����ɋʂ�����
	Vector3d bl_x = Vector3d(0, 0, 50e-3 + 1.001e-3);
	Vector3d bl_v = Vector3d(0, 0, 0);
	Vector3d bl_w = Vector3d(0, 0, 1);		// -�� �����Ɏ��]
	Quaterniond q = Quaterniond(1, 0, 0, 0); // �����N�H�[�^�j�I���͑S�� [0,0,0,1] �Ƃ���D�i�������R���X�g���N�^�̎d�l�㏇�Ԃ��قȂ�j

	this->BL->set_y(bl_x, bl_v, q, bl_w);
	Vector3d Fbi, Tbi, Fib, Tib, Fs, Ts;

	// �׏d�v�Z���C�X���C�X�̌��ʂ��擾
	this->BNP.get_FT(Fbi, Tbi, Tib, Fs, Ts);
	this->BNP.save(FO.BNP[0]);

	// 1. �e�X���C�X�Ђ̂��ׂ薀�C�̐������m�F
	// �ڐG�ȉ~�� -�� ���ɂ���X���C�X�Ђ̂��ׂ薀�C�� -�� ����
	for (int i = 0; i < 11; i++) {
		double fs_xai = FO.BNP[0].GV[0].SL[i].fs_[0];
		double ps_zeta = FO.BNP[0].GV[0].SL[i].ps_[2];
		EXPECT_LT(fs_xai, 0);
	}
	// �ڐG�ȉ~�� +�� ���ɂ���X���C�X�Ђ̂��ׂ薀�C�� +�� ����
	for (int i = 12; i < 21; i++) {
		double fs_xai = FO.BNP[0].GV[0].SL[i].fs_[0];
		double ps_zeta = FO.BNP[0].GV[0].SL[i].ps_[2];
		EXPECT_GT(fs_xai, 0);
	}
	// 2. ���C�� << �w���c�׏d �~ ���C�W�� �ɂȂ��Ă��邩�m�F
	Vector3d Fn = -Fbi - Fs;
	EXPECT_LT(Fs.norm(), Fn.norm() * 0.1);

	// 3. �X���C�X10 �̉e�������������΁C�ʎ��]�Ƃ͋t�����ɖ��C���������邱�Ƃ��m�F
	Vector3d Fs__, Ts__, ts_, fs_;
	// �X���C�X10�̖��C�͗������W�n�ŏo�͂����̂ŁC�������W�n�ɒ����Ă���o��
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


// (�É��step0) calc_force_case1 �Ɠ��������ŉ׏d�v�Z���C���ʂ��m�F
TEST_F(BS_BallCylinderPairTest, get_F0_case1) {
	this->BS_init();
	// �ʐڋߗ� 1um�C�ʂ̐i�s�����͗����ڐ������ɂȂ�悤�ɐݒ�
	Vector3d bl_x = Vector3d(0, 0, 50e-3 + 1.001e-3);
	Vector3d bl_v = Vector3d(0.5, 0.8660254, 0);
	Vector3d bl_w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // �����N�H�[�^�j�I���͑S�� [0,0,0,1] �Ƃ���D�i�������R���X�g���N�^�̎d�l�㏇�Ԃ��قȂ�j

	this->BL->set_y(bl_x, bl_v, q, bl_w);

	Vector2d Fbc;
	Vector3d Fcb, Tcb;
	this->BNP.get_F0(Fbc, Fcb, Tcb);

	// 1. �ʂ���i�b�g�ɓ����w���c�׏d�� -�� �����i= +z �����j�̐��������ɂȂ��Ă��邩�m�F
	double err = 1e-3;
	EXPECT_NEAR(Fbc[0], Fbc.norm(), Fbc.norm() * err);
	EXPECT_NEAR(Fcb[2], Fcb.norm(), Fcb.norm() * err);
	return;
}

// (�É��step2) calc_force_case1 �Ɠ��������ŉ׏d�v�Z���C���ʂ��m�F
TEST_F(BS_BallCylinderPairTest, get_F1_case1) {
	this->BS_init();
	// �ʐڋߗ� 1um�C�ʂ̐i�s�����͗����ڐ������ɂȂ�悤�ɐݒ�
	Vector3d bl_x = Vector3d(0, 0, 50e-3 + 1.001e-3);
	Vector3d bl_v = Vector3d(0.5, 0.8660254, 0);
	Vector3d bl_w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // �����N�H�[�^�j�I���͑S�� [0,0,0,1] �Ƃ���D�i�������R���X�g���N�^�̎d�l�㏇�Ԃ��قȂ�j

	this->BL->set_y(bl_x, bl_v, q, bl_w);

	Vector2d Fbc;
	Vector3d Fcb, Tcb;
	this->BNP.get_F1(true, Fbc, Fcb, Tcb);

	// �X�e�b�v1 �̐ڐG�׏d�̓X�e�b�v0�̉׏d�ɖ��C�׏d�����������̂Ȃ̂ŁC
	// �X�e�b�v0�̌��ʂɗ\�z����門�C�׏d�����Z���ăX�e�b�v1�̗\���l�𓱏o�D
	Vector3d Fcb_ex = Vector3d(-3.1495739376147817 * cos30, -3.1495739376147817 * sin30, 31.495739376147817);
	EXPECT_LT((Fcb - Fcb_ex).norm(), 1e-2);
	return;
}
// (�É��step2) �V���t�g���ŉ׏d�v�Z���C���ʂ��m�F
TEST_F(BS_BallCylinderPairTest, get_F1_case2) {
	this->BS_init();
	// �ʐڋߗ� 1um�C�ʂ̐i�s�����͗����ڐ������ɂȂ�悤�ɐݒ�
	Vector3d bl_x = Vector3d(0, 0, 50e-3 - 1.001e-3);
	Vector3d bl_v = Vector3d(0.5, 0.8660254, 0);
	Vector3d bl_w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // �����N�H�[�^�j�I���͑S�� [0,0,0,1] �Ƃ���D�i�������R���X�g���N�^�̎d�l�㏇�Ԃ��قȂ�j

	this->BL->set_y(bl_x, bl_v, q, bl_w);
	Vector2d Fbc;
	Vector3d Fcb, Tcb;
	this->BSP.get_F1(false, Fbc, Fcb, Tcb);

	// �X�e�b�v1 �̐ڐG�׏d�̓X�e�b�v0�̉׏d�ɖ��C�׏d�����������̂Ȃ̂ŁC
	// �X�e�b�v0�̌��ʂɗ\�z����門�C�׏d�����Z���ăX�e�b�v1�̗\���l�𓱏o�D
	Vector3d Fcb_ex = Vector3d(-3.0082504011194946 * cos30, -3.0082504011194946 * sin30, -30.082504011194946);
	EXPECT_LT((Fcb - Fcb_ex).norm(), 1e-2);
	return;
}

// (�É��step3) calc_force_case1 �Ɠ��������ŉ׏d�v�Z���C���ʂ��m�F
TEST_F(BS_BallCylinderPairTest, get_F2_case1) {
	this->BS_init();
	// �ʐڋߗ� 1um�C�ʂ̐i�s�����͗����ڐ������ɂȂ�悤�ɐݒ�
	Vector3d bl_x = Vector3d(0, 0, 50e-3 + 1.001e-3);
	Vector3d bl_v = Vector3d(sin30, cos30, 0);
	Vector3d bl_w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // �����N�H�[�^�j�I���͑S�� [0,0,0,1] �Ƃ���D�i�������R���X�g���N�^�̎d�l�㏇�Ԃ��قȂ�j

	this->BL->set_y(bl_x, bl_v, q, bl_w);
	Vector2d Fbc;
	Vector3d Fcb, Tcb;
	this->BNP.get_F1(true, Fbc, Fcb, Tcb);

	Vector3d vF, vT;
	this->BNP.get_F2(vF, vT);

	// �X�e�b�v2 �̐ڐG�׏d�̓X�e�b�v0�̉׏d�ɂ��ׂ葬�x�����������̂ɂȂ��Ă��邩�m�F
	Vector3d vF_ex = Vector3d(31.495739376147817 * sin30, 31.495739376147817 * cos30, 0);
	EXPECT_LT((vF - vF_ex).norm(), 1e-2);
	return;
}