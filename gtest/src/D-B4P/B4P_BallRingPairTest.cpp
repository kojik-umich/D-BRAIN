#include "pch.h"
#include "B4P_StabIn.h"
#include "B4P_StabOut.h"

class B4P_BallRingPairTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	B4P_BallInnerRingPair BIP;
	B4P_BallOuterRingPair BOP;
	Ball BL;
	B4P_InnerRing IR;
	B4P_OuterRing OR;
	B4P_StabIn FI;
	B4P_StabOut FO;

	// ���󏏌���50BSWZ02�̂��̂ɐݒ�
	void init_50BSWZ() {
		// ����
		FI.LB.eta0 = 0.5;	// �S�x�iPa*s�j(�K��)
		FI.LB.beta0 = 0.5;	// ���x�S�x�W���i�K���j
		FI.LB.k0 = 0.145;	// ���M�`����[W/(m�EK)]
		FI.LB.alpha0 = 0.2;	// ���͔S�x�W��[mm2/kgf]
		FI.LB.lm0 = -1;	// ���j�X�J�X����


		// �ʕ����l
		double D = 0.0111;				// �ʌa[m]
		double E = 208000000000;		// �����O��[Pa]
		double por = 0.29;				// �|�A�\����[-]
		double den = 7830;				// ���x[kg/m^3]
		double rms = 0.00000002;		// �e��rms[m]
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}

		this->BL.init(D, E, por, den, rms, x_const, Rx_const);
		int msmax = 21;
		FO.allocate(1, _MAX_CONTACT_, msmax);

		// ���֕����l
		FI.ballpcd = 0.0675;
		double clr = 0.000008;
		FI.IR.R[0] = 0.00589;
		FI.IR.R[1] = 0.00589;
		FI.IR.Rox[0] = 0.000167;
		FI.IR.Rox[1] = -0.000167;
		// �aR���SPCD[m]
		FI.IR.Rod[0] = FI.ballpcd - clr * 0.5 + 2 * sqrt((FI.IR.R[0] - D * 0.5)*(FI.IR.R[0] - D * 0.5) - FI.IR.Rox[0] * FI.IR.Rox[0]);
		FI.IR.Rod[1] = FI.ballpcd - clr * 0.5 + 2 * sqrt((FI.IR.R[1] - D * 0.5)*(FI.IR.R[1] - D * 0.5) - FI.IR.Rox[1] * FI.IR.Rox[1]);

		FI.IR.hedge[0];
		FI.IR.hedge[1];
		FI.IR.rms = 0.00000006;
		FI.IR.m = 1;
		FI.IR.E = 208000000000;
		FI.IR.por = 0.29;
		FI.IR.Ix, FI.IR.Iyz, FI.IR.Iyz;
		FI.BIP.mu = 0.1;
		FI.BIP.dzeta = 0.2;
		// �������_
		FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::RollingResistanceNothing;
		FI.TB.coulomb = B4P_In::Tribology::Tangent;
		FI.TB.coulomb_slope = 1e100;
		FI.TB.hysteresis = B4P_In::Tribology::HysteresisNothing;
		this->IR.init(FI.IR, FI.ballpcd);
		this->BIP.link(&this->BL, &this->IR);
		this->BIP.init(FI);

		// �O�֕����l
		FI.OR.Rox[0] = 0.0001665;
		FI.OR.Rox[1] = -0.0001665;
		FI.OR.R[0] = 0.00589;
		FI.OR.R[1] = 0.00589;
		// �aR���SPCD[m]
		FI.OR.Rod[0] = FI.ballpcd + clr * 0.5 - 2 * sqrt((FI.OR.R[0] - D * 0.5)*(FI.OR.R[0] - D * 0.5) - FI.OR.Rox[0] * FI.OR.Rox[0]);
		FI.OR.Rod[1] = FI.ballpcd + clr * 0.5 - 2 * sqrt((FI.OR.R[1] - D * 0.5)*(FI.OR.R[1] - D * 0.5) - FI.OR.Rox[1] * FI.OR.Rox[1]);
		FI.OR.hedge[0];
		FI.OR.hedge[1];
		FI.OR.rms = 0.00000006;
		FI.OR.m = 1;
		FI.OR.E = 207760000000;
		FI.OR.por = 0.29;
		FI.OR.Ix, FI.OR.Iyz, FI.OR.Iyz;
		FI.LB.eta0 = 0.5;	// �S�x�iPa*s�j(�K��)
		FI.LB.beta0 = 0.5;	// ���x�S�x�W���i�K���j
		FI.LB.k0 = 0.145;	// ���M�`����[W/(m�EK)]
		FI.LB.alpha0 = 0.2;	// ���͔S�x�W��[mm2/kgf]
		FI.LB.lm0 = -1;	// ���j�X�J�X����
		FI.BOP.mu = 0.1;	// ���C�W��
		FI.BOP.dzeta = 0.2;	// �����W��
		this->OR.init(FI.OR, FI.ballpcd);
		this->BOP.link(&this->BL, &this->OR);
		this->BOP.init(FI);


		return;
	};

	// ���󏏌���25BSWZ01�̂��̂ɐݒ�
	void init_25BSWZ() {
		// ���̓f�[�^�͉��L�̂��̂�p����
		// \\ans00978\kiken\BRAIN\�\�[�X�R�[�h��\DBRAINpp\docs\06_D4Bv100\05_�P�̃e�X�g\B4P_BallRingPair
		// ����^�ԁ@25BSWZ01
		// Pair�N���X�����l
		this->FI.BOP.dzeta = 0.2;	// ������
		this->FI.BOP.mu = 0.1;		// ���C�W��
		this->FI.BIP.dzeta = 0.2;	// ������
		this->FI.BIP.mu = 0.1;		// ���C�W��

		// �ʕ����l
		double D = 0.00635;				// �ʌa[m]
		double E = 207760000000;		// �����O��[Pa]
		double por = 0.29;				// �|�A�\����[-]
		double den = 7830;				// ���x[kg/m^3]
		double rms = 0.00000002;		// �e��rms[m]
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}
		this->BL.init(D, E, por, den, rms, x_const, Rx_const);
		int msmax = 21;
		FO.allocate(1, _MAX_CONTACT_, msmax);

		// ���֕����l
		FI.ballpcd = 0.0355;
		double clr = 0;
		FI.IR.Rox[0] = 0.0000725;
		FI.IR.Rox[1] = -0.0000725;
		FI.IR.R[0] = 0.00332105;
		FI.IR.R[1] = 0.00332105;
		// �aR���SPCD[m]�i�apcd���a�C0.0357535 m ���炢�j
		FI.IR.Rod[0] = FI.ballpcd - clr * 0.5 + 2 * sqrt((FI.IR.R[0] - D * 0.5)*(FI.IR.R[0] - D * 0.5) - FI.IR.Rox[0] * FI.IR.Rox[0]);
		FI.IR.Rod[1] = FI.ballpcd - clr * 0.5 + 2 * sqrt((FI.IR.R[1] - D * 0.5)*(FI.IR.R[1] - D * 0.5) - FI.IR.Rox[1] * FI.IR.Rox[1]);
		FI.IR.hedge[0];
		FI.IR.hedge[1];
		FI.IR.rms = 0.00000006;
		FI.IR.m = 1;
		FI.IR.E = 207760000000;
		FI.IR.por = 0.29;
		FI.IR.Ix, FI.IR.Iyz, FI.IR.Iyz;

		// �������_
		FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::RollingResistanceNothing;
		FI.TB.coulomb = B4P_In::Tribology::Tangent;
		FI.TB.coulomb_slope = 1000;
		FI.TB.hysteresis = B4P_In::Tribology::HysteresisNothing;
		this->IR.init(FI.IR, FI.ballpcd);
		this->BIP.link(&this->BL, &this->IR);
		this->BIP.init(FI);

		// �O�֕����l
		FI.OR.Rox[0] = 0.0000725;
		FI.OR.Rox[1] = -0.0000725;
		FI.OR.R[0] = 0.00332105;
		FI.OR.R[1] = 0.00332105;
		FI.OR.Rod[0] = FI.ballpcd + clr * 0.5 - 2 * sqrt((FI.OR.R[0] - D * 0.5)*(FI.OR.R[0] - D * 0.5) - FI.OR.Rox[0] * FI.OR.Rox[0]);
		FI.OR.Rod[1] = FI.ballpcd + clr * 0.5 - 2 * sqrt((FI.OR.R[1] - D * 0.5)*(FI.OR.R[1] - D * 0.5) - FI.OR.Rox[1] * FI.OR.Rox[1]);
		FI.OR.hedge[0];
		FI.OR.hedge[1];
		FI.OR.rms = 0.00000006;
		FI.OR.m = 1;
		FI.OR.E = 207760000000;
		FI.OR.por = 0.29;
		FI.OR.Ix, FI.OR.Iyz, FI.OR.Iyz;
		FI.LB.eta0 = 0.5;	// �S�x�iPa*s�j(�K��)
		FI.LB.beta0 = 0.5;	// ���x�S�x�W���i�K���j
		FI.LB.k0 = 0.145;	// ���M�`����[W/(m�EK)]
		FI.LB.alpha0 = 0.2;	// ���͔S�x�W��[mm2/kgf]
		FI.LB.lm0 = -1;	// ���j�X�J�X����
		this->OR.init(FI.OR, FI.ballpcd);
		this->BOP.link(&this->BL, &this->OR);
		this->BOP.init(FI);

	};

	virtual void SetUp() {
		Vector3d x = Vector3d(0, 0, 0);
		Vector3d v = Vector3d(0, 0, 0);
		Quaterniond q = Quaterniond(1, 0, 0, 0); // �����N�H�[�^�j�I���͑S�� [0,0,0,1] �Ƃ���D�i�������R���X�g���N�^�̎d�l�㏇�Ԃ��قȂ�j
		Vector3d w = Vector3d(0, 0, 0);
		FI.msmax = 21;
		Rigid::l = 1;
		Rigid::t = 1;
		Rigid::g = Vector3d(0, 0, 0);
		this->OR.set_y(x, v, q, w);
		this->IR.set_y(x, v, q, w);
		int msmax = 21;
		FO.allocate(1, _MAX_CONTACT_, msmax);
	}

	virtual void TearDown() {

	}
};



// (1) �ێ���E�ʁE���O�֎��]���x���p�����[�^�Ƃ��������쐬���C��r�D
TEST_F(B4P_BallRingPairTest, get_us_ur_1) {
	this->init_25BSWZ();
	double wc = 44.4604328912874;		// �ێ�����]���x[rad/s]
	double wrn_0 = 0.0;					// �O�։�]���x[rad/s]
	double wrn_1 = 104.719755119660;	// ���։�]���x[rad/s]
	double omega_x = -223.906546549487;	// �ʎ��]���x ��x[rad/s]
	double omega_z = 287.884790174018;  // �ʎ��]���x ��z[rad/s]
	double dx_out = 4.389357066845456E-006;		// �O�֑��ڋߗ�(balac�v�Z����)
	double dx_in = 4.518845377961239E-006;		// ���֑��ڋߗ�(balac�v�Z����)
	double sin_alp1 = 0.537033675917149;		// �O�֐ڐG�p�̐���(balac�v�Z����)
	double sin_alp2 = 0.537125019737040;		// ���֐ڐG�p�̐���(balac�v�Z����)
	double cos_alp1 = 0.843560804525029;		// �O�֐ڐG�p�̗]��(balac�v�Z����)
	double cos_alp2 = -0.843502645622694;		// ���֐ڐG�p�̗]��(balac�v�Z����)

	// �ʑ��Έʒu�͐ڋߗʂ���t�Z
	double R_out = this->OR.GV[0].r - this->BL.r;
	double R_in = this->IR.GV[0].r - this->BL.r;
	Vector3d er = Vector3d(0, 0, 1);			// �a�����x�N�g��
	Vector3d bl_x_in = Vector3d(0, 0, 0);
	Vector3d bl_x_out = Vector3d(0, 0, 0);

	Vector3d x = Vector3d(0, 0, 0);
	Vector3d v = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // �����N�H�[�^�j�I���͑S�� [0,0,0,1] �Ƃ���D�i�������R���X�g���N�^�̎d�l�㏇�Ԃ��قȂ�j
	Vector3d wo = Vector3d(wrn_0, 0, 0);
	this->OR.set_y(x, v, q, wo);
	Vector3d wi = Vector3d(wrn_1, 0, 0);
	this->IR.set_y(x, v, q, wi);
	Vector3d bl_x = Vector3d(0, 0, 0.5 * FI.ballpcd);
	Vector3d bl_v = Vector3d(0, -0.5 * FI.ballpcd * wc, 0);
	Vector3d bl_w = Vector3d(omega_x, 0, omega_z);
	this->BL.set_y(bl_x, bl_v, q, bl_w);

	// �ʒ��S�ɑ΂���ڐG�_�̑��΋���
	Vector3d _p_out = Vector3d(this->BL.r * sin_alp1, 0, this->BL.r * cos_alp1 + 0.5 * FI.ballpcd);
	Vector3d _p_in = Vector3d(this->BL.r * sin_alp2, 0, this->BL.r * cos_alp2 + 0.5 * FI.ballpcd);
	Vector3d uso_, uro_, usi_, uri_;
	this->BOP.get_us_ur(_p_out, uso_, uro_);

	double uso_norm = uso_.norm(), uro_norm = uro_.norm();
	this->BIP.get_us_ur(_p_in, usi_, uri_);
	double usi_norm = usi_.norm(), uri_norm = uri_.norm();

	// baltac�Ɏ�������Ă��������Q�l�Ɏ��삵����(�ꕔ�ϐ��̒�`���@���قȂ邽�߉���)
	// �i�ڐG�_�̑��x�̌v�Z���@���{�v���O�����̌v�Z���@�ɍ��킹����ł���j
	double u1, u2, baltac_uso, baltac_usi, baltac_Uo = 0, baltac_Ui = 0;
	u1 = 0.5 * FI.ballpcd  * wc - 0.5 * (FI.ballpcd + this->BL.D * cos_alp1) * wrn_0;
	u2 = 0.5 * this->BL.D *(-omega_x * cos_alp1 + omega_z * sin_alp1);
	baltac_uso = u1 - u2;
	baltac_Uo = (u1 + u2) * 0.5;
	u1 = 0.5 * FI.ballpcd  * wc - 0.5 * (FI.ballpcd + this->BL.D * cos_alp2) * wrn_1;
	u2 = 0.5 * this->BL.D *(-omega_x * cos_alp2 + omega_z * sin_alp2);
	baltac_usi = u1 - u2;
	baltac_Ui = (u1 + u2) * 0.5;

	double err = 0.05;

	cout << uro_norm << "\t" << abs(baltac_Uo) << endl;
	cout << uso_norm << "\t" << abs(baltac_uso) << endl;
	cout << uri_norm << "\t" << abs(baltac_Ui) << endl;
	cout << usi_norm << "\t" << abs(baltac_usi) << endl;
	// ���e�덷
	EXPECT_NEAR(uro_norm, abs(baltac_Uo), abs(baltac_Uo)*err);
	EXPECT_NEAR(uso_norm, abs(baltac_uso), abs(baltac_uso)*err);
	EXPECT_NEAR(uri_norm, abs(baltac_Ui), abs(baltac_Ui)*err);
	EXPECT_NEAR(usi_norm, abs(baltac_usi), abs(baltac_usi)*err);

	return;
}

// (1) �ێ���E�ʁE���O�֎��]���x���p�����[�^�Ƃ��������쐬���C��r�D
TEST_F(B4P_BallRingPairTest, get_us_ur2_1) {

	this->init_25BSWZ();
	double wc = 44.4604328912874;		// �ێ�����]���x[rad/s]
	double wrn_0 = 0.0;					// �O�։�]���x[rad/s]
	double wrn_1 = 104.719755119660;	// ���։�]���x[rad/s]

	double omega_x = -223.906546549487;	// �ʎ��]���x ��x[rad/s]
	double omega_z = 287.884790174018;  // �ʎ��]���x ��z[rad/s]
	double dx_out = 4.389357066845456E-006;		// �O�֑��ڋߗ�(balac�v�Z����)
	double dx_in = 4.518845377961239E-006;		// ���֑��ڋߗ�(balac�v�Z����)
	double sin_alp1 = 0.537033675917149;		// �O�֐ڐG�p�̐���(balac�v�Z����)
	double sin_alp2 = 0.537125019737040;		// ���֐ڐG�p�̐���(balac�v�Z����)
	double cos_alp1 = 0.843560804525029;		// �O�֐ڐG�p�̗]��(balac�v�Z����)
	double cos_alp2 = -0.843502645622694;		// ���֐ڐG�p�̗]��(balac�v�Z����)

	// �ʑ��Έʒu�͐ڋߗʂ���t�Z
	double R_out = this->OR.GV[0].r - this->BL.r;
	double R_in = this->IR.GV[0].r - this->BL.r;
	Vector3d er = Vector3d(0, 0, 1);			// �a�����x�N�g��
	Vector3d bl_x_in = Vector3d(0, 0, 0.5 * FI.ballpcd);
	Vector3d bl_x_out = Vector3d(0, 0, 0.5 * FI.ballpcd);

	Vector3d x = Vector3d(0, 0, 0);
	Vector3d v = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // �����N�H�[�^�j�I���͑S�� [0,0,0,1] �Ƃ���D�i�������R���X�g���N�^�̎d�l�㏇�Ԃ��قȂ�j
	Vector3d wo = Vector3d(wrn_0, 0, 0);
	this->OR.set_y(x, v, q, wo);
	Vector3d wi = Vector3d(wrn_1, 0, 0);
	this->IR.set_y(x, v, q, wi);
	Vector3d bl_v = Vector3d(0, -0.5 * FI.ballpcd * wc, 0);
	Vector3d bl_w = Vector3d(omega_x, 0, omega_z);
	this->BL.set_y(bl_x_in, bl_v, q, bl_w);

	// �ʒ��S�ɑ΂���ڐG�_�̑��΋���
	Vector3d _p_out = Vector3d(this->BL.r * sin_alp1, 0, this->BL.r * cos_alp1 + 0.5 * FI.ballpcd);
	Vector3d _p_in = Vector3d(this->BL.r * sin_alp2, 0, this->BL.r * cos_alp2 + 0.5 * FI.ballpcd);
	Vector3d uso_, uro_, usi_, uri_;
	this->BOP.get_us_ur2(_p_out, er, uso_, uro_);

	double uso_norm = uso_.norm(), uro_norm = uro_.norm();
	this->BIP.get_us_ur2(_p_in, er, usi_, uri_);
	double usi_norm = usi_.norm(), uri_norm = uri_.norm();

	// baltac�Ɏ�������Ă��������Q�l�Ɏ��삵����(�ꕔ�ϐ��̒�`���@���قȂ邽�߉���)
	// �i�ڐG�_�̑��x�̌v�Z���@���{�v���O�����̌v�Z���@�ɍ��킹����ł���j
	double u1, u2, baltac_uso, baltac_usi, baltac_Uo = 0, baltac_Ui = 0;
	u1 = 0.5 * (FI.ballpcd + this->BL.D * cos_alp1)*(wc - wrn_0);
	u2 = 0.5 * this->BL.D *(-(omega_x - (wc - wrn_0))* cos_alp1 + omega_z * sin_alp1);
	baltac_uso = u1 - u2;
	baltac_Uo = (u1 + u2) * 0.5;
	u1 = 0.5 * (FI.ballpcd + this->BL.D * cos_alp2)*(wc - wrn_1);
	u2 = 0.5 * this->BL.D *(-(omega_x - (wc))* cos_alp2 + omega_z * sin_alp2);
	baltac_usi = u1 - u2;
	baltac_Ui = (u1 + u2) * 0.5;

	double err = 0.05;

	cout << uro_norm << "\t" << abs(baltac_Uo) << endl;
	cout << uso_norm << "\t" << abs(baltac_uso) << endl;
	cout << uri_norm << "\t" << abs(baltac_Ui) << endl;
	cout << usi_norm << "\t" << abs(baltac_usi) << endl;
	// ���e�덷
	EXPECT_NEAR(uro_norm, abs(baltac_Uo), abs(baltac_Uo)*err);
	EXPECT_NEAR(uso_norm, abs(baltac_uso), abs(baltac_uso)*err);
	EXPECT_NEAR(uri_norm, abs(baltac_Ui), abs(baltac_Ui)*err);
	EXPECT_NEAR(usi_norm, abs(baltac_usi), abs(baltac_usi)*err);

	return;
}


// (1) ���C�W�����Œ�l�ɂ����Ƃ��ɋʁ[���O�֊Ԃ̉׏d����v�Z�Ɠ����ɂȂ��Ă��邩����
TEST_F(B4P_BallRingPairTest, calc_force_1) {
	this->init_25BSWZ();

	// �ڐG���
	double sin_alp1 = 0.540517404112816;		// �O�֐ڐG�p�̐���(balac�v�Z����)
	double cos_alp1 = 0.841332832980588;		// �O�֐ڐG�p�̗]��(balac�v�Z����)
	double dx_out = 4.389357066845456E-006;		// �O�֑��ڋߗ�(balac�v�Z����)
	double _Qo = 21.3811947318042;				// �]���̉׏d[kgf](balac�v�Z����)

	this->BOP.TR = new Tribology::Stab_Traction();	// �g���N�V�����W����0.2�ɌŒ�
	this->BOP.CL = new Tribology::Stab_Coulomb();	// �N�[�������C�W����0.3�ɌŒ�
	this->BOP.RR = new Tribology::Stab_RollingResist();	// �]����S����R��萔�ɌŒ�
	this->BOP.FT = new Tribology::FilmThicknessNothing();	// ����������0�ɌŒ�

	// �ʂ�����]�ŕ��i�^�����Ă���ꍇ������
	int i = 1;			// �ڐG����a�ԍ�
	Vector3d x = Vector3d(0, 0, 0), v = Vector3d(0, 0, 0), w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	this->OR.set_y(x, v, q, w);
	Vector3d er = Vector3d(0, 0, 1);
	double dx = this->OR.GV[i].Rx;
	double R_out = this->OR.GV[i].r + dx_out - this->BL.r;			// �aR���S����ʒ��S�̋���
	double rp = this->OR.GV[i].r + dx_out * 0.5;					// �aR���S����ڐG�_�̋���
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + dx, 0, cos_alp1 * R_out) + er * this->OR.GV[i].Rr;
	Vector3d bl_v = Vector3d(0, -1, 0);
	this->BL.set_y(bl_x_out, bl_v, q, w);


	Vector3d Fbi, Tbi, Fib, Tib;

	this->BOP.calc_force(Fbi, Tbi, Fib, Tib);

	// ��v�Z�̌��ʂƔ�r
	double _fn = _Qo * 9.8;
	Vector3d _Fn = Vector3d(-sin_alp1, 0, -cos_alp1) * _fn;	// �����׏d[N]
	Vector3d _Fs = Vector3d(0, _fn * 0.3, 0);				// ���薀�C[N]
	Vector3d _Fbi = _Fn + _Fs;
	double err = 0.05;
	EXPECT_NEAR((Fbi - _Fbi).norm(), 0, _Fbi.norm() * err);
	double _rb = this->BL.r - dx_out * 0.5;
	double _Ts_norm = _fn * 0.3 * _rb;
	Vector3d _Ts = Vector3d(-cos_alp1, 0, sin_alp1) * _Ts_norm;			// ���薀�C[N]
	Vector3d _Tr = Vector3d(cos_alp1, 0, -sin_alp1) * 0.1;
	Vector3d _Tbi = _Ts + _Tr;
	EXPECT_NEAR((Tbi - _Tbi).norm(), 0, _Tbi.norm() * err);

	// write()�̃e�X�g�����˂ĐڐG�ȉ~��ڐG�p������

	this->BOP.write(FO.BOP[0]);
	Vector3d Fn(FO.BOP[0].GV[i].Fn);		// �����׏d[N]
	EXPECT_NEAR((Fn - _Fn).norm(), 0, _Fn.norm() * err);
	Vector3d Fs(FO.BOP[0].GV[i].Fs);		// ���薀�C[N]
	EXPECT_NEAR((Fs - _Fs).norm(), 0, _Fs.norm() * err);
	double ea = FO.BOP[0].GV[i].a;											// �ȉ~�����a[mm](����)
	double _ea = 0.668610936573505e-3;										// �ȉ~�����a[mm](���Ғl)
	EXPECT_NEAR(_ea, ea, _ea * err);
	double eb = FO.BOP[0].GV[i].b;											// �ȉ~�Z���a[mm](����)
	double _eb = 9.720461846754701e-5;								// �ȉ~�Z���a[mm](���Ғl)
	EXPECT_NEAR(_eb, eb, _eb * err);
	double alp = FO.BOP[0].GV[i].phi;											// �ڐG�p[��](����)
	double _alp = Unit::deg2rad(32.7188677338554);										// �ڐG�p[��](���Ғl)
	EXPECT_NEAR(_alp, alp, _alp * err);

	return;
}

// (1) �ڐG�p0���ɂ�����ʐڋߗʂƊe������x�N�g�����m�F
TEST_F(B4P_BallRingPairTest, how_Contact_1) {
	this->init_25BSWZ();

	// �ڐG�p0���E�ڋߗ�1mm�̂Ƃ��̋ʒ��S�ʒu���t�Z
	// �ʒ��Sy���W = �ڋߗ� + �aR���a - �ʔ��a + �aR���Sy���W
	// �������C�a���S�ʒu��(0.0000725, 0.035246430660370774* 0.5, 0)
	// �ڐG����a�� +X ���Ƃ���

	Vector3d x(-0.0000725, (0.00332105 + 0.035246430660370774 * 0.5 - 0.00635 * 0.5 + 0.001), 0);

	Vector3d er, eg;
	double dx;

	int i = 1; 
	bool c = this->BOP.how_Contact(i, x, er, eg, dx);

	double err = 0.01;					// ���e�덷
	Vector3d expcted_er(0, 1, 0);		// �ʈʑ��x�N�g��
	Vector3d expcted_eg(0, 1, 0);		// �ʐڐG�����x�N�g��
	double expected_dx = 0.001;			// �ʐڋߗ�

	EXPECT_NEAR((er - expcted_er).norm(), 0, er.norm()*err);
	EXPECT_NEAR((eg - expcted_eg).norm(), 0, eg.norm()*err);
	EXPECT_NEAR(dx, expected_dx, dx*err);

	return;
};

// (2) �{�[�����S���g�[���X����o�Ă�����ڐG�Ȃ��̔���ƂȂ邱�Ƃ��m�F
TEST_F(B4P_BallRingPairTest, how_Contact_2) {
	this->init_25BSWZ();

	// �{�[�����S���g�[���X���O���ɐݒ�C�ڐG����a��+X���ɐݒ�
	Vector3d x(-0.0000725, 1.0, 0);
	int i = 1;

	Vector3d er, eg;
	double dx;
	bool c = this->BOP.how_Contact(i, x, er, eg, dx);

	EXPECT_EQ(false, c);
	return;
};

// (3) �H�����ݗʂ����l�̏ꍇ�ڐG�Ȃ��̔���ƂȂ邱�Ƃ��m�F
TEST_F(B4P_BallRingPairTest, how_Contact_3) {
	this->init_25BSWZ();

	// �ڋߗ�-10um���ڐG�p0���ƂȂ�悤�ɋʈʒu���w��
	Vector3d x(-0.0000725, (0.00332105 + 0.035246430660370774 * 0.5 - 0.00635 * 0.5 - 0.0001), 0);
	//Vector3d v(0.0, -0.01, 0.0);
	int i = 1;

	Vector3d er, eg;
	double dx;

	bool c = this->BOP.how_Contact(i, x,  er, eg, dx);
	EXPECT_EQ(false, c);
	return;
};


// (4) �ʂ��a���S�_�Əd�Ȃ�Ƃ��C�ڐG�Ȃ��ɂȂ邩�m�F
TEST_F(B4P_BallRingPairTest, how_Contact_4) {
	this->init_25BSWZ();

	// �ʈʒu���a���S�ʒu�ɔz�u
	Vector3d x(-0.0000725, 0.035246430660370774 * 0.5, 0);
	int i = 1;

	Vector3d er, eg;
	double dx;

	// �ʏ�͂��肦�Ȃ����C�ʌa��1.2�{�ɂ��ċʂ��a���S�ɂ��鎞�ł��O�ւƐڐG����悤�ɒ���
	double D = 0.00635 * 1.2;		// �ʌa[m]
	double E = 207760000000;		// �����O��[Pa]
	double por = 0.29;				// �|�A�\����[-]
	double den = 7830;				// ���x[kg/m^3]
	double rms = 0.00000002;		// �e��rms[m]
	bool x_const[3], Rx_const[3];
	for (int i = 0; i < 3; i++) {
		x_const[i] = false;
		Rx_const[i] = false;
	}
	this->BL.init(D, E, por, den, rms, x_const, Rx_const);


	bool c = this->BOP.how_Contact(i, x, er, eg, dx);

	double err = 0.05; // ���e�덷
	Vector3d expcted_er(0, 1, 0);
	Vector3d expcted_eg(0, 1, 0);
	double expected_dx = 0.001;

	EXPECT_EQ(false, c);
	return;
};

// (5) �ʑ��p�ƐڐG�p����[�����C�ڐG�E��ڐG�̋��E�t�߂̋������m�F
TEST_F(B4P_BallRingPairTest, how_Contact_5) {
	this->init_25BSWZ();

	// �ʂ��ʑ��p60���C�ڐG�p30���C�ڋߗ�0�̈ʒu�ɔz�u
	// x���W = �aR���Sx���W + (�ڋߗ� + �aR���a - �ʔ��a) * sin30��
	// y���W = ((�ڋߗ� + �aR���a - �ʔ��a) * cos30�� + �aPCD���a) * sin60��
	// z���W = ((�ڋߗ� + �aR���a - �ʔ��a) * cos30�� + �aPCD���a) * cos60��

	Vector3d x(5.25E-07, 0.01537169, 0.008874849);
	int i = 1;
	Vector3d er, eg;
	double dx;

	// �ڋߗ�0�̈ʒu���� +x, +y, +z �����ɔ����ʈړ�������ƐڐG���邩����
	Vector3d ep_x(1e-6, 0, 0), ep_y(0, 1e-6, 0), ep_z(0, 0, 1e-6);
	bool cx1 = this->BOP.how_Contact(i, x + ep_x, er, eg, dx);
	bool cy1 = this->BOP.how_Contact(i, x + ep_y, er, eg, dx);
	bool cz1 = this->BOP.how_Contact(i, x + ep_z, er, eg, dx);
	EXPECT_EQ(true, cx1);
	EXPECT_EQ(true, cy1);
	EXPECT_EQ(true, cz1);

	// �ڋߗ�0�̈ʒu���� -x, -y, -z �����ɔ����ʈړ�������Ɣ�ڐG�ɂȂ邩����
	bool cx2 = this->BOP.how_Contact(i, x - ep_x, er, eg, dx);
	bool cy2 = this->BOP.how_Contact(i, x - ep_y, er, eg, dx);
	bool cz2 = this->BOP.how_Contact(i, x - ep_z, er, eg, dx);
	EXPECT_EQ(false, cx2);
	EXPECT_EQ(false, cy2);
	EXPECT_EQ(false, cz2);
	return;
};


// (1) �ʑ��x��ς��Čv�Z���s���C�������l�����Ȃ��ꍇ�Ƃ̔�r���s���D
TEST_F(B4P_BallRingPairTest, calc_DynamicHertz_1) {
	this->init_25BSWZ();


	/*****d4b���͏���(���󏏌���SetUp�֐��Őݒ�)*****/
	// �ʂ����������猩�� +y �������ɂ���ꍇ��z��D
	// �ʑ��x�� +y �����i�O�֊O���ɐi�s�j

	int i = 0;									// �a�ԍ�
	double sin_alp1 = 0.540517404112816;		// �O�֐ڐG�p�̐���(balac�v�Z����)
	double cos_alp1 = 0.841332832980588;		// �O�֐ڐG�p�̗]��(balac�v�Z����)
	Vector3d er= Vector3d(0, 1, 0);					// �a�����x�N�g��
	Vector3d eg_out= Vector3d(sin_alp1, cos_alp1, 0);	// �a���S���{�[�������x�N�g��(baltac�ڐG�p����t�Z)
	double dx_out= 4.389357066845456E-006;		// �O�֑��ڋߗ�(balac�v�Z����)
	double Rx_out, Ry_out, cos_alp, sin_alp, a_out, b_out;
	Vector3d p_out;
	// �ʑ��Έʒu�͐ڋߗʂ���t�Z
	// �ʂ��ʑ��p 90�� (+y����)�C�ڐG�p 30���C�ڋߗ� 4.389 um �̈ʒu�ɔz�u
	// x���W = �aR���Sx���W + (�ڋߗ� + �aR���a - �ʔ��a) * sin30��
	// y���W = ((�ڋߗ� + �aR���a - �ʔ��a) * cos30�� + �aPCD���a) * sin90��
	// z���W = ((�ڋߗ� + �aR���a - �ʔ��a) * cos30�� + �aPCD���a) * cos90��
	double R_out = this->OR.GV[0].r + dx_out - this->BL.r;
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + this->OR.GV[0].Rx, cos_alp1 * R_out, 0) + er * this->OR.GV[0].Rr;
	
	// ��r�ΏۂƂȂ錸���Ȃ��̏ꍇ���Ɍv�Z�D
	double k_out, k_in;
	this->BOP.how_Contact(i, bl_x_out, er, eg_out, dx_out);
	this->BOP.calc_Hertz(i, bl_x_out, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out, k_out);
	double _Qout = k_out * pow(dx_out, 1.5);
	// (1) �ʂ��ڐG�ʂɋ߂Â������ɐi�s����ꍇ�C�ڐG�͂��傫���Ȃ邩�m�F
	Vector3d bl_v1 = Vector3d(1, 0, 0);		// �ʑ��x�i�ǖʂ���߂Â������j
	double Qout1 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v1, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_GE(Qout1, _Qout);	

	// (2) �ʂ��ڐG�ʂ��痣�������ɐi�s����ꍇ�C�ڐG�͂��������Ȃ邩�m�F
	Vector3d bl_v2 = Vector3d(-1, 0, 0);	// �ʑ��x�i�ǖʂ��痣�������j
	double Qout2 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v2, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_LE(Qout2, _Qout);

	// (3) �ʂ��ڐG�ʂƕ��s�Ȍ����ɐi�s����ꍇ�C�ڐG�͂��ς��Ȃ����Ƃ��m�F(z�����̂�)
	Vector3d bl_v3 = Vector3d(0, 0, -1);	// �ʑ��x�i�ǖʕ��s�����j
	double Qout3 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v3, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_NEAR(Qout3, _Qout, 1e-6);

	// (4) �ʂ��ڐG�ʂƕ��s�Ȍ����ɐi�s����ꍇ�C�ڐG�͂��ς��Ȃ����Ƃ��m�F(xy�����̂�)
	Vector3d bl_v4 = Vector3d(0.841332832980588, -0.540517404112816, 0);	// �ʑ��x�i�ǖʕ��s�����j
	double Qout4 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v4, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_NEAR(Qout4, _Qout, 1e-6);

	// (5) �����͂�Hertz�ڐG�͂̍��v�l�����̎��C0�ɕ␳���Ă��邱�Ƃ��m�F
	Vector3d bl_v5 = Vector3d(-1e9, 0, 0);		// �ʑ��x�i�ǖʂ��痣�������j
	double Qout5 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v5, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_NEAR(Qout5, 0, 1e-6);	

	return;
};




// (1) �����������Ȃ��Ƃ��̃w���c�ڐG�͂��v�Z���Cbaltac�̌v�Z���ʂƔ�r
// baltac�Ɠ����̐ڋߗʁE�ڐG�p����͂��������ŁC�]���̉׏d�E�ȉ~���a��]��
// �i���͏����̓A�L�V�A���׏d�̏ꍇ��p����D�����͍l�����Ȃ��D�j
TEST_F(B4P_BallRingPairTest, calc_Hertz_1) {
	this->init_25BSWZ();

	// ���󏏌���SetUp�֐��Őݒ�
	// d4b���͏���
	int i = 0;									// �a�ԍ�
	double sin_alp1 = 0.540517404112816;		// �O�֐ڐG�p�̐���(balac�v�Z����)
	double sin_alp2 = 0.540517404112816;		// ���֐ڐG�p�̐���(balac�v�Z����)
	double cos_alp1 = 0.841332832980588;		// �O�֐ڐG�p�̗]��(balac�v�Z����)
	double cos_alp2 = -0.841332832980588;		// ���֐ڐG�p�̗]��(balac�v�Z����)
	Vector3d er = Vector3d(0, 1, 0);			// �a�����x�N�g��
	Vector3d eg_out = Vector3d(sin_alp1, cos_alp1, 0); // �a���S���{�[�������x�N�g��(baltac�ڐG�p����t�Z)
	Vector3d eg_in = Vector3d(sin_alp2, cos_alp2, 0); // �a���S���{�[�������x�N�g��(baltac�ڐG�p����t�Z)
	double dx_out = 4.389357066845456E-006;		// �O�֑��ڋߗ�(balac�v�Z����)
	double dx_in = 4.518845377961239E-006;		// ���֑��ڋߗ�(balac�v�Z����)
	double Rx_out, Ry_out, Rx_in, Ry_in, cos_alp, sin_alp, a_out, b_out, k_out, a_in, b_in, k_in;
	Vector3d p_out, p_in;
	// �ʑ��Έʒu�͐ڋߗʂ���t�Z
	double R_out = this->OR.GV[0].r + dx_out - this->BL.r;
	double R_in = this->IR.GV[0].r + dx_in - this->BL.r;
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + this->OR.GV[0].Rx, cos_alp1 * R_out, 0) + er * this->OR.GV[0].Rr;
	Vector3d bl_x_in = Vector3d(sin_alp2 * R_in + this->IR.GV[0].Rx, cos_alp2 * R_in, 0) + er * this->IR.GV[0].Rr;
	// �e�X�g�Ώۊ֐��̌Ăяo��
	this->BOP.calc_Hertz(i, bl_x_out, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out, k_out);
	double Qout = k_out * pow(dx_out, 1.5);
	this->BIP.calc_Hertz(i, bl_x_in, er, eg_in, dx_in, Rx_in, Ry_in, p_in, cos_alp, sin_alp, a_in, b_in, k_in);
	double Qin = k_in * pow(dx_in, 1.5);

	// ��r�ΏۂƂȂ�balac�v�Z���ʁi�f�o�b�O���[�h�ɂĎ擾�j
	// 1: �O�ցC2:����
	double ea1 = 0.668610936573505;				// �ȉ~�����a[mm]
	double ea2 = 0.685074685553727;				// �ȉ~�����a[mm]
	double eb1 = 9.720461846754701E-002;		// �ȉ~�Z���a[mm]
	double eb2 = 8.272063984522840E-002;		// �ȉ~�Z���a[mm]
	double _Q1 = 21.3811947318042;				// �]���̉׏d[kgf]
	double _Q2 = 21.3811947318041;				// �]���̉׏d[kgf]
	double _rho_0 = 2 / 6.35;					// �]���̋ȗ�[1/mm]
	double _rho3_in = -0.301109588834856;		// �a�ȗ�[1/mm]
	double _rho2_in = 5.579585936574204E-002;	// �O���̋ȗ�(���]�O���ڐG���a�̋t���j[1/mm]
	double _Rx_in = 1 / (_rho_0 + _rho2_in);
	double _Ry_in = 1 / (_rho_0 + _rho3_in);

	double _rho3_out = -0.301109588834856;
	double _rho2_out = -4.119892685701449E-002;
	double _Rx_out = 1 / (_rho_0 + _rho2_out);
	double _Ry_out = 1 / (_rho_0 + _rho3_out);
	// �ڐG�ʒu(�e�X�g�P�[�X)�͐ڋߗʂƐڐG�p����v�Z
	Vector3d _p_out = Vector3d(sin_alp1 * (this->OR.GV[0].r + dx_out * 0.5) + this->OR.GV[0].Rx, cos_alp1 * (this->OR.GV[0].r + dx_out * 0.5), 0)
		+ er * this->OR.GV[0].Rr;
	Vector3d _p_in = Vector3d(sin_alp2 * (this->IR.GV[0].r + dx_in * 0.5) + this->IR.GV[0].Rx, cos_alp2 * (this->IR.GV[0].r + dx_in * 0.5), 0)
		+ er * this->IR.GV[0].Rr;


	cout << "���O�֑�����" << endl;
	cout << "�]���̉׏d�F" << Qout << "\tbaltac�F" << _Q1 * 9.8 << endl;
	cout << "�ڐG�ȉ~�����a�F" << a_out << "\tbaltac�F" << ea1 / 1000 << endl;
	cout << "�ڐG�ȉ~�Z���a�F" << b_out << "\tbaltac�F" << eb1 / 1000 << endl;
	cout << "�ȗ�Rx�F" << Rx_out << "\tbaltac�F" << _Rx_out / 1000 << endl;
	cout << "�ȗ�Ry�F" << Ry_out << "\tbaltac�F" << _Ry_out / 1000 << endl;
	cout << "�ڐG�_�ʒu�덷�F" << (p_out - _p_out).norm() << endl;
	//cout << "p = " << p_out << "\thand�F" << _p_out <<endl;
	cout << "�����֑�����" << endl;
	cout << "�]���̉׏d�F" << Qin << "\tbaltac�F" << _Q2 * 9.8 << endl;
	cout << "�ڐG�ȉ~�����a�F" << a_in << "\tbaltac�F" << ea2 / 1000 << endl;
	cout << "�ڐG�ȉ~�Z���a�F" << b_in << "\tbaltac�F" << eb2 / 1000 << endl;
	cout << "�ȗ�Rx�F" << Rx_in << "\tbaltac�F" << _Rx_in / 1000 << endl;
	cout << "�ȗ�Ry�F" << Ry_in << "\tbaltac�F" << _Ry_in / 1000 << endl;
	cout << "�ڐG�_�ʒu�덷�F" << (p_in - _p_in).norm() << endl;
	//cout << "p = " << p_in << "\thand�F" << _p_in << endl;
	// baltac�ł�Brew-Hamrock�̎��ł͂Ȃ��C
	// ���������W�J�ɂ��ߎ��ŋ��߂Ă��邽�ߊ��S�Ɉ�v���Ȃ��D
	// �ibaltac������K(k')��E(k')�̌v�Z�͉^���͊w�Up12�Ɍf�ځj
	// �덷��2%�ȓ��ł���ΉƂ����D
	double err = 0.02;
	EXPECT_NEAR(Qout, _Q1*9.8, Qout*err);
	EXPECT_NEAR(Qin, _Q2*9.8, Qin*err);
	EXPECT_NEAR(a_out, ea1 / 1000, a_out*err);
	EXPECT_NEAR(a_in, ea2 / 1000, a_in*err);
	EXPECT_NEAR(b_out, eb1 / 1000, b_out*err);
	EXPECT_NEAR(b_in, eb2 / 1000, b_in*err);
	EXPECT_NEAR(Rx_out, _Rx_out / 1000, Rx_out*err);
	EXPECT_NEAR(Rx_in, _Rx_in / 1000, Rx_in*err);
	EXPECT_NEAR(Ry_out, _Ry_out / 1000, Ry_out*err);
	EXPECT_NEAR(Ry_in, _Ry_in / 1000, Ry_in*err);
	EXPECT_NEAR((p_out - _p_out).norm(), 0, p_out.norm()*err);
	EXPECT_NEAR((p_in - _p_in).norm(), 0, p_in.norm()*err);
	return;
};

// (1)�ʂ����O�։~�������Ɉړ������Ƃ��̊��薀�C���v�Z���C��v�Z�̐��l�Ɣ�r
// �������C�g���N�V�����W���E�N�[�������C�W���͒萔�l�ɌŒ�
TEST_F(B4P_BallRingPairTest, calc_Sliding_1) {
	this->init_25BSWZ();

	int i = 0;	// �a�ԍ�
	double a = 0.668610936573505e-3; // �ڐG�ȉ~���a[m]
	double Pmax = 1;		// �ő�ʈ��i�g��Ȃ��j
	double F_norm = 100;	// �]���̉׏d[N]
	double fratio = 0.4;	// �����ڐG����
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double sin_alp1 = 0.540517404112816;		// �O�֐ڐG�p�̐���(balac�v�Z����)
	double sin_alp2 = 0.540517404112816;		// ���֐ڐG�p�̐���(balac�v�Z����)
	double cos_alp1 = 0.841332832980588;		// �O�֐ڐG�p�̗]��(balac�v�Z����)
	double cos_alp2 = -0.841332832980588;		// ���֐ڐG�p�̗]��(balac�v�Z����)
	double dx_out = 4.389357066845456E-006;		// �O�֑��ڋߗ�(balac�v�Z����)


	this->BOP.TR = new Tribology::Stab_Traction();	// �g���N�V�����W����0.2�ɌŒ�
	this->BOP.CL = new Tribology::Stab_Coulomb();	// �N�[�������C�W����0.3�ɌŒ�

	// �ʂ�Y�����i�������j�ɖ���]�Ŋ����Ă����Ԃ�����
	Vector3d x = Vector3d(0, 0, 0);
	Vector3d v = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d w = Vector3d(0, 0, 0);
	this->OR.set_y(x, v, q, w);
	Vector3d bl_v = Vector3d(0, 1, 0);
	Vector3d er = Vector3d(0, 0, 1);
	double R_out = this->OR.GV[0].r + dx_out - this->BL.r;
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + this->OR.GV[0].Rx, 0, cos_alp1 * R_out) + er * this->OR.GV[0].Rr;
	Vector3d p = Vector3d(sin_alp1 * (this->OR.GV[0].r + dx_out * 0.5)+ this->OR.GV[0].Rx, 0, cos_alp1 * (this->OR.GV[0].r + dx_out * 0.5)) + er * this->OR.GV[0].Rr;
	this->BL.set_y(bl_x_out, bl_v, q, w);

	this->BOP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);

	// ��v�Z�̌��ʂƔ�r
	double _fs = 100 * (0.4 * 0.2 + 0.6 * 0.3);
	double bp = this->BL.r - dx_out * 0.5;
	Vector3d _Fbs = Vector3d(0, -_fs, 0);
	Vector3d _Tbs = Vector3d(_fs *bp * cos_alp1, 0, -_fs * bp * sin_alp1);
	double err = 0.05;
	EXPECT_NEAR((Fbs - _Fbs).norm(), 0, _Fbs.norm()*err);
	EXPECT_NEAR((Tbs - _Tbs).norm(), 0, _Tbs.norm()*err);


	return;
}

// (2) �ʂ��O�ւɑ΂��ăX�s�������ɉ�]�������̊��薀�C���v�Z
// �������C�N�[�������C�W���݂̂Ƃ���
TEST_F(B4P_BallRingPairTest, calc_Sliding_2) {
	this->init_25BSWZ();

	int i = 0;	// �a�ԍ�
	double a = 0.668610936573505e-3; // �ڐG�ȉ~���a[m]
	double Pmax = 1;		// �ő�ʈ��i�g��Ȃ��j
	double F_norm = 100;	// �]���̉׏d[N]
	double fratio = 0.0;	// �����ڐG����
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double sin_alp1 = 0.540517404112816;		// �O�֐ڐG�p�̐���(balac�v�Z����)
	double cos_alp1 = 0.841332832980588;		// �O�֐ڐG�p�̗]��(balac�v�Z����)
	double dx_out = 4.389357066845456E-006;		// �O�֑��ڋߗ�(balac�v�Z����)


	this->BOP.TR = new Tribology::Stab_Traction();	// �g���N�V�����W����0.2�ɌŒ�
	this->BOP.CL = new Tribology::Stab_Coulomb();	// �N�[�������C�W����0.3�ɌŒ�

	// �ʂ����x0�ŁC�ڐG�ʐ������������ɉ�]���Ă���ꍇ������
	Vector3d x = Vector3d(0, 0, 0);
	Vector3d v = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d w = Vector3d(0, 0, 0);
	this->OR.set_y(x, v, q, w);
	Vector3d er = Vector3d(0, 0, 1);
	double R_out = this->OR.GV[0].r + dx_out - this->BL.r;
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + this->OR.GV[0].Rx, 0, cos_alp1 * R_out) + er * this->OR.GV[0].Rr;
	Vector3d bl_w = Vector3d(sin_alp1, 0, cos_alp1);
	Vector3d p = Vector3d(sin_alp1 * (this->OR.GV[0].r + dx_out * 0.5) + this->OR.GV[0].Rx, 0, cos_alp1 * (this->OR.GV[0].r + dx_out * 0.5)) + er * this->OR.GV[0].Rr;
	this->BL.set_y(bl_x_out, v, q, bl_w);
	this->BOP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);



	// ��v�Z�̌��ʂƔ�r
	double _fs = 100 * (0.4 * 0.2 + 0.6 * 0.3);
	double bp = this->BL.r - dx_out * 0.5;
	Vector3d _Fbs = Vector3d(0, 0, 0);
	EXPECT_NEAR((Fbs - _Fbs).norm(), 0, 0.1);			// ���薀�C�̃g�[�^����0
	Vector3d eg = Vector3d(sin_alp1, 0, cos_alp1);
	double dir = Tbs.dot(eg);
	Vector3d Dir = Tbs.cross(eg);
	EXPECT_LT(dir, 0);								// �g���N�̌������ڐG�ʂɑ΂��ă{�[�����ɂȂ��Ă��邱�Ƃ��m�F
	EXPECT_NEAR(Dir.norm(), 0, 0.1);				// �g���N���ڐG�ʂɑ΂��Đ����ł��邩�m�F
	return;
}

// (3) �y�O��+���z�ʐڐG�p��30���łقڏ��]���肵�Ă���Ƃ��̖��C�͂̌�����]��
TEST_F(B4P_BallRingPairTest, calc_Sliding_3) {
	this->init_50BSWZ();
	int i = 1;							// �a�ԍ�
	double a = 0.3160272198e-3;			// �ڐG�ȉ~���a[m]
	double Pmax = 2e6;					// �ő�ʈ��i�g��Ȃ��j
	double F_norm = 20;		// �]���̉׏d[N]
	double fratio = 0.0;				// �����ڐG����
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double sin_alp1 = 0.5;				// 30����sin
	double cos_alp1 = sqrt(3) * 0.5;	// 30����cos
	double dx_out = 0;		// �O�֑��ڋߗ�(balac�v�Z����)

	// ���̏����F�Ȃ��C�N�[�������C�Ftan�J�[�u
	this->BOP.TR = new Tribology::Stab_Traction();
	this->BOP.CL = new Tribology::Tangent();
	this->BOP.FT = new Tribology::FilmThicknessNothing();
	// �ʂ����]����(���v���)�ŁC�ڐG�ʐ������������ɉ�]���Ă���ꍇ������
	double wx = 20;
	Vector3d bl_v = Vector3d(0, wx * (-dx_out + this->BL.r), 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d bl_w = Vector3d(wx * cos_alp1, 0, -wx * sin_alp1);

	// �ʈʒu�͍aR���S�̈ʒu����Ɍ���
	double R_out = this->OR.GV[i].r + dx_out - this->BL.r;
	Vector3d er = Vector3d(0, 0, 1);
	Vector3d Ro = Vector3d(this->OR.GV[i].Rx, 0, this->OR.GV[i].Rr);
	Vector3d bl_x = Ro + Vector3d(sin_alp1 * R_out, 0, cos_alp1 * R_out);
	this->BL.set_y(bl_x, bl_v, q, bl_w);
	Vector3d zero = Vector3d::Zero();
	this->OR.set_y(zero, zero, q, zero);
	double rr = this->OR.GV[i].r + dx_out * 0.5;
	Vector3d p = Ro + Vector3d(sin_alp1 * rr, 0, cos_alp1 * rr);

	this->BOP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);



	// ���C�́E���[�����g�̕]��
	Vector3d e_fs = Fbs.normalized();
	Vector3d e_ts = Tbs.normalized();
	Vector3d e_fsex(0, -1, 0);
	Vector3d e_tsex(cos_alp1, 0, -sin_alp1);
	EXPECT_NEAR((e_fs - e_fsex).norm(), 0, 1e-6);			// ���薀�C�̌���
	EXPECT_NEAR((e_ts - e_tsex).norm(), 0, 1e-6);			// ���[�����g�̌���

	return;
}


// (4) �y�O��+���z�ʐڐG�p��0���ŋʎ��]����x�����-45���������Ă���ꍇ�̋ʃ��[�����g�̌������������v�Z�ł��Ă��邩����
TEST_F(B4P_BallRingPairTest, calc_Sliding_4) {
	this->init_50BSWZ();
	int i = 1;							// �a�ԍ�
	double a = 0.3160272198e-3;			// �ڐG�ȉ~���a[m]
	double Pmax = 2e6;					// �ő�ʈ��i�g��Ȃ��j
	double F_norm = 20;					// �]���̉׏d[N]
	double fratio = 0.0;				// �����ڐG����
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double dx_out = 0;					// �O�֑��ڋߗ�
	double cos_45 = 1.0 / sqrt(2.0); 
	double sin_45 = 1.0 / sqrt(2.0);

	this->BOP.TR = new Tribology::Stab_Traction();
	this->BOP.CL = new Tribology::Tangent();
	this->BOP.FT = new Tribology::FilmThicknessNothing();
	// �ʂ����x0�ŁC�ڐG�ʐ������������ɉ�]���Ă���ꍇ������
	double wx = 100;
	Vector3d bl_v = Vector3d(0, wx * (this->BL.r + dx_out * 0.5) * cos_45, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d bl_w = Vector3d(wx * sin_45, 0,  -wx * cos_45);

	// �ʈʒu�͍aR���S�̈ʒu����Ɍ���
	double R_out = this->OR.GV[i].r + dx_out - this->BL.r;
	Vector3d bl_x = Vector3d(this->OR.GV[i].Rx, 0, R_out + this->OR.GV[i].Rr);
	this->BL.set_y(bl_x, bl_v, q, bl_w);
	Vector3d zero = Vector3d::Zero();
	this->OR.set_y(zero, zero, q, zero);
	double rr = this->OR.GV[i].r + dx_out * 0.5;
	Vector3d p = Vector3d(this->OR.GV[i].Rx, 0, rr + this->OR.GV[i].Rr);


	this->BOP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);

	// ��v�Z�̌��ʂƔ�r
	EXPECT_GT(Tbs.z(), 0);								// z+�����Ƀ��[�����g������
	EXPECT_LT(Tis.z(), 0);								// z-�����Ƀ��[�����g������

	return;
}

// (5) �y�O��+���z�ʐڐG�p��0���ŋʎ��]����x�����+45���������Ă���ꍇ�̋ʃ��[�����g�̌������������v�Z�ł��Ă��邩����
TEST_F(B4P_BallRingPairTest, calc_Sliding_5) {
	this->init_50BSWZ();
	int i = 1;							// �a�ԍ�
	double a = 0.3160272198e-3;			// �ڐG�ȉ~���a[m]
	double Pmax = 2e6;					// �ő�ʈ��i�g��Ȃ��j
	double F_norm = 20;		// �]���̉׏d[N]
	double fratio = 0.0;				// �����ڐG����
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double dx_out = 0;		// �O�֑��ڋߗ�
	double cos_45 = 1.0 / sqrt(2.0);
	double sin_45 = 1.0 / sqrt(2.0);

	this->BOP.TR = new Tribology::Stab_Traction();
	this->BOP.CL = new Tribology::Tangent();
	this->BOP.FT = new Tribology::FilmThicknessNothing();
	// �ʂ����x0�ŁC�ڐG�ʐ������������ɉ�]���Ă���ꍇ������
	double wx = 100;
	Vector3d bl_v = Vector3d(0, wx * (this->BL.r + dx_out * 0.5) * cos_45, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d bl_w = Vector3d(wx * sin_45, 0, wx * cos_45);

	// �ʈʒu�͍aR���S�̈ʒu����Ɍ���
	double R_out = this->OR.GV[i].r + dx_out - this->BL.r;
	Vector3d bl_x = Vector3d(this->OR.GV[i].Rx, 0, R_out + this->OR.GV[i].Rr);
	this->BL.set_y(bl_x, bl_v, q, bl_w);
	Vector3d zero = Vector3d::Zero();
	this->OR.set_y(zero, zero, q, zero);
	double rr = this->OR.GV[i].r + dx_out * 0.5;
	Vector3d p = Vector3d(this->OR.GV[i].Rx, 0, rr + this->OR.GV[i].Rr);


	this->BOP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);
	// ��v�Z�̌��ʂƔ�r
	EXPECT_LT(Tbs.z(), 0);								// z-�����Ƀ��[�����g������

	return;
}

// (6) �y����-���z�ʐڐG�p��0���ŋʎ��]����x�����+45���������Ă���ꍇ�̋ʃ��[�����g�̌������������v�Z�ł��Ă��邩����
TEST_F(B4P_BallRingPairTest, calc_Sliding_6) {
	this->init_50BSWZ();
	int i = 0;							// �a�ԍ�
	double a = 0.3160272198e-3;			// �ڐG�ȉ~���a[m]
	double Pmax = 2e6;					// �ő�ʈ��i�g��Ȃ��j
	double F_norm = 20;		// �]���̉׏d[N]
	double fratio = 0.0;				// �����ڐG����
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double dxi = 0;		// ���֑��ڋߗ�
	double cos_45 = 1.0 / sqrt(2.0);
	double sin_45 = 1.0 / sqrt(2.0);
	double cos_m135 = -1.0 / sqrt(2.0);
	double sin_m135 = -1.0 / sqrt(2.0);
	this->BIP.TR = new Tribology::Stab_Traction();
	this->BIP.CL = new Tribology::Tangent();
	this->BIP.FT = new Tribology::FilmThicknessNothing();
	// �ڐG�_�ɂ����銊�葬�x��0�̎�������
	double wx = 100;
	Vector3d bl_v = Vector3d(0, wx * (this->BL.r + dxi * 0.5) * cos_45, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d bl_w = Vector3d(wx * sin_45, 0, wx * cos_45);

	// �ʈʒu�͍aR���S�̈ʒu����Ɍ���
	double R_in = this->IR.GV[i].r + dxi - this->BL.r;
	Vector3d bl_x = Vector3d(this->IR.GV[i].Rx, 0, this->IR.GV[i].Rr - R_in);
	this->BL.set_y(bl_x, bl_v, q, bl_w);
	Vector3d zero = Vector3d::Zero();
	this->IR.set_y(zero, zero, q, zero);
	double rr = this->IR.GV[i].r + dxi * 0.5;
	Vector3d p = Vector3d(this->IR.GV[i].Rx, 0, rr + this->IR.GV[i].Rr);
	this->BIP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);
	EXPECT_LT(Tbs.z(), 0);								// z-�����Ƀ��[�����g������

	return;
}


// (1) ���̐ڐG�p�̎��̊O�֐ڐG�p���v�Z���D��v�Z�̐��l�Ɣ�r
TEST_F(B4P_BallRingPairTest, ContactAngle_1) {
	this->init_25BSWZ();

	double err = 0.01;
	double sin_alp1 = 0.5;				// �O�֐ڐG�p�̐���
	double cos_alp1 = 0.8660254;		// �O�֐ڐG�p�̗]��
	double alp = this->BOP.ContactAngle(cos_alp1, sin_alp1);
	double _alp = Unit::deg2rad(30);
	EXPECT_NEAR(_alp, alp, alp * err);
	return;
}

// (2) ���̐ڐG�p�̎��̊O�֐ڐG�p���v�Z���D��v�Z�̐��l�Ɣ�r
TEST_F(B4P_BallRingPairTest, ContactAngle_2) {
	this->init_25BSWZ();

	double err = 0.01;
	double sin_alp1 = -0.5;				// �O�֐ڐG�p�̐���
	double cos_alp1 = 0.8660254;		// �O�֐ڐG�p�̗]��
	double alp = this->BOP.ContactAngle(cos_alp1, sin_alp1);
	double _alp = Unit::deg2rad(-30);
	EXPECT_NEAR(_alp, alp, abs(alp) * err);
	return;
}

// (3) ���̐ڐG�p�̎��̓��֐ڐG�p���v�Z���D��v�Z�̐��l�Ɣ�r
TEST_F(B4P_BallRingPairTest, ContactAngle_3) {
	this->init_25BSWZ();

	double err = 0.01;
	double sin_alp1 = 0.5;				// �O�֐ڐG�p�̐���
	double cos_alp1 = -0.8660254;		// �O�֐ڐG�p�̗]��
	double alp = this->BIP.ContactAngle(cos_alp1, sin_alp1);
	double _alp = Unit::deg2rad(30);
	EXPECT_NEAR(_alp, alp, alp * err);
	return;
}

// (4) ���̐ڐG�p�̎��̓��֐ڐG�p���v�Z���D��v�Z�̐��l�Ɣ�r
TEST_F(B4P_BallRingPairTest, ContactAngle_4) {
	this->init_25BSWZ();

	double err = 0.01;
	double sin_alp1 = -0.5;				// �O�֐ڐG�p�̐���
	double cos_alp1 = -0.8660254;		// �O�֐ڐG�p�̗]��
	double alp = this->BIP.ContactAngle(cos_alp1, sin_alp1);
	double _alp = Unit::deg2rad(-30);
	EXPECT_NEAR(_alp, alp, abs(alp) * err);
	double alp__ = Unit::rad2deg(atan2(sin_alp1, cos_alp1));
	return;
}

// (1) ��������̊m�F
TEST_F(B4P_BallRingPairTest, init_Tribology_1) {
	this->init_25BSWZ();

	FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::Aihara;
	this->BOP.init_Tribology(FI.TB);
	Tribology::AiharaR* aih = dynamic_cast<Tribology::AiharaR*>(this->BOP.RR);
	EXPECT_TRUE(aih != nullptr);

	FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::Fujiwara;
	this->BOP.init_Tribology(FI.TB);
	Tribology::Fujiwara* fj = dynamic_cast<Tribology::Fujiwara*>(this->BOP.RR);
	EXPECT_TRUE(fj != nullptr);

	FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::Houpert;
	this->BOP.init_Tribology(FI.TB);
	Tribology::Houpert* hou = dynamic_cast<Tribology::Houpert*>(this->BOP.RR);
	EXPECT_TRUE(hou != nullptr);

	FI.TB.coulomb = B4P_In::Tribology::Coulomb::CoulombNothing;
	this->BOP.init_Tribology(FI.TB);
	Tribology::CoulombNothing* cln = dynamic_cast<Tribology::CoulombNothing*>(this->BOP.CL);
	EXPECT_TRUE(cln != nullptr);

	FI.TB.coulomb = B4P_In::Tribology::Coulomb::Tangent;
	this->BOP.init_Tribology(FI.TB);
	Tribology::Tangent* tan = dynamic_cast<Tribology::Tangent*>(this->BOP.CL);
	EXPECT_TRUE(tan != nullptr);

	FI.TB.filmThickness = B4P_In::Tribology::FilmThickness::FilmThicknessNothing;
	this->BOP.init_Tribology(FI.TB);
	Tribology::FilmThicknessNothing* flt = dynamic_cast<Tribology::FilmThicknessNothing*>(this->BOP.FT);
	EXPECT_TRUE(flt != nullptr);

	FI.TB.filmThickness = B4P_In::Tribology::FilmThickness::HamrockDowsonHc;
	this->BOP.init_Tribology(FI.TB);
	Tribology::HamrockDowsonHc* hmd = dynamic_cast<Tribology::HamrockDowsonHc*>(this->BOP.FT);
	EXPECT_TRUE(hmd != nullptr);
	return;
}

