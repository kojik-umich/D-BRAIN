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

		// ���֕����l
		this->D = 0.00635;				// �ʌa[m]
		FI.ballpcd = 0.0355;
		double clr = 0;
		FI.IR.Rox[0] = 0.0000725;
		FI.IR.Rox[1] = -0.0000725;
		FI.IR.R[0] = 0.00332105;
		FI.IR.R[1] = 0.00332105;
		// �aR_pcd (= 0.035753569339629220)
		FI.IR.Rod[0] = FI.ballpcd + (FI.IR.R[0] - D * 0.5) * 0.866025403784439;
		FI.IR.Rod[1] = FI.ballpcd + (FI.IR.R[1] - D * 0.5) * 0.866025403784439;
		FI.IR.hedge[0];
		FI.IR.hedge[1];
		FI.IR.rms = 0.00000006;
		FI.IR.m = 1;
		FI.IR.E = 207760000000;
		FI.IR.por = 0.29;
		FI.IR.Ix, FI.IR.Iyz, FI.IR.Iyz;

		// �������_
		FI.TB.rollingresistance;
		this->IR.init(FI.IR, FI.ballpcd);

		// �O�֕����l
		FI.OR.Rox[0] = 0.0000725;
		FI.OR.Rox[1] = -0.0000725;
		FI.OR.R[0] = 0.00332105;
		FI.OR.R[1] = 0.00332105;
		// �aR_pcd (= 0.035246430660370774)
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
		Quaterniond q = Quaterniond(1, 0, 0, 0); // �����N�H�[�^�j�I���͑S�� [0,0,0,1] �Ƃ���D�i�������R���X�g���N�^�̎d�l�㏇�Ԃ��قȂ�j
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


// (1) �e�X���C�X�Ђ̒��S�ƍaR���S��e������x�N�g���̑��Ί֌W��������������
TEST_F(B4P_OuterRingTest, calc_slice_test1) {

	Vector3d p(this->OR.GV[0].Rx, 0, this->OR.GV[0].r + this->OR.GV[0].Rr);
	double a = 1e-4;
	int ms = 21;
	Vector3d ps[21];

	this->OR.calc_slice(p, a, ms, 0, this->D, ps);
	EXPECT_NEAR((ps[20] - ps[0]).norm(), 2 * a, 1e-1 * a) << "�ڐG�ȉ~�̒[����[�̋������C�ڐG�ȉ~�̒����ɓ��������m�F";

	// �ڐG�ȉ~�̒[����[�ƁC�aR���S����ڐG�_�ւ̃x�N�g�������s���Ă��邩�̂��̊m�F
	Vector3d OP(0, 0, this->OR.GV[0].r);
	double ip1 = OP.dot((ps[20] - ps[0]));
	EXPECT_NEAR(ip1, 0 , 1e-6);

	// �ڐG�ȉ~�̒[����[�ƁC�������x�N�g�����������Ă��邩�̊m�F
	Vector3d xai(0, 1, 0);
	double ip2 = xai.dot(ps[20] - ps[0]);
	EXPECT_NEAR(ip2, 0, 1e-6);

	// �X���C�X�В��S����aR���S�ւ̋��������ꂼ�ꓙ�������̊m�F
	Vector3d O(this->OR.GV[0].Rx, 0, this->OR.GV[0].Rr);
	double l1 = (ps[20] - O).norm();
	double l2 = (ps[0] - O).norm();
	EXPECT_NEAR(l1, l2, 1e-6);

	//  �ڐG�ȉ~�̕В[����aR���S�ւ̋��� l1���ڐG�_����aR���S�ւ̋��� l0 ���傫�����̊m�F
	// (�d�l�����̂��ߕۗ�)
	double l0 = (ps[10] - O).norm();
	EXPECT_GE(l0, l1);

	// �ڐG�ȉ~�̕В[����a���S�ւ̋������C�a���a��1%�ȓ��̍��ɂƂǂ܂��Ă��邩�̊m�F
	Vector3d Gc(this->OR.GV[0].Rx, 0, this->OR.GV[0].Rr);
	double OP0 = (ps[20] - Gc).norm();
	double r = this->OR.GV[0].r;
	EXPECT_NEAR(l1, r, 1e-5) << "�ڐG�ȉ~�̕В[����a���S�ւ̋������C�a���a��1%�ȓ��̍��ɂƂǂ܂��Ă��邩�̊m�F";
	return;
};


// (2) �񎟊֐��ߎ��ł͂Ȃ��~�ʂ̎���p���Čv�Z�C��r
TEST_F(B4P_OuterRingTest, calc_slice_test2) {
	// �ڐG�p30���̎��̐ڐG�_
	Vector3d Ro(this->OR.GV[0].Rx, 0, this->OR.GV[0].Rr);	// ���ڂ���aR���S
	double theta = Numeric::pi * 0.166666;					// �ڐG�p
	Vector3d p = Ro + this->OR.GV[0].r * Vector3d(sin(theta), 0, cos(theta));

	// �v���O�����̎��Ōv�Z
	double a = 1e-6;					// �ڐG�ȉ~[m]
	int ms = 21;						// ������
	Vector3d *ps = new Vector3d[ms];	// �ڐG�_���W[m]
	this->OR.calc_slice(p, a, ms, 0, this->D, ps);


	// �~�̎���p���Čv�Z
	Vector3d *ps_ex = new Vector3d[ms];	// �ڐG�_���W�i���Ғl�j[m]
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

// (3) ���ւƊO�ւŋʂ�ڐG�p30���ɂȂ�悤�ɋ��񂾎��C���ւƊO�ւ̃X���C�X�Ђ����s�ɂȂ邩���؁D
TEST_F(B4P_OuterRingTest, calc_slice_test3) {
	// ����x���W = 2 * (�aR���a - �ʔ��a) * sin30�� - Rx * 2
	double bl_r = 0.00635 * 0.5;
	double Ox = 2 * (this->IR.GV[0].r - bl_r) * 0.5 - (this->IR.GV[0].Rx) * 2;
	Vector3d x_out	= Vector3d(0, 0, 0);
	Vector3d x_in	= Vector3d(Ox, 0, 0);
	Vector3d v		= Vector3d(0, 0, 0);
	Quaterniond q	= Quaterniond(1, 0, 0, 0); // �����N�H�[�^�j�I���͑S�� [0,0,0,1] �Ƃ���D�i�������R���X�g���N�^�̎d�l�㏇�Ԃ��قȂ�j
	Vector3d w		= Vector3d(0, 0, 0);

	this->OR.set_y(x_out, v, q, w);
	this->IR.set_y(x_in, v, q, w);

	double a = 1e-6;								// �ڐG�ȉ~[m]
	int ms = 21;									// ������
	Vector3d *ps_out = new Vector3d[ms];			// �ڐG�_���W[m]
	Vector3d p_out = Vector3d(this->OR.GV[1].Rx, 0, this->OR.GV[1].Rr) + this->OR.GV[1].r * Vector3d(0.5, 0, 0.866025403784439);
	this->OR.calc_slice(p_out, a, ms, 1, this->D, ps_out);
	
	Vector3d p_in_other = Vector3d(Ox + this->IR.GV[0].Rx, 0, this->IR.GV[0].Rr) + this->IR.GV[0].r * Vector3d(-0.5, 0, -0.866025403784439);
	Vector3d *ps_in = new Vector3d[ms];				// �ڐG�_���W[m]
	Vector3d p_in  = p_out + 2 * bl_r * Vector3d(-0.5, 0, -0.866025403784439);
	this->IR.calc_slice(p_in_other, a, ms, 0, this->D, ps_in);
	Vector3d pdist_out = ps_out[20] - ps_out[0];
	Vector3d pdist_in = ps_in[20] - ps_in[0];

	Vector3d cross = pdist_out.normalized().cross(pdist_in.normalized());
	EXPECT_NEAR(cross.norm(), 0, 1e-3);
	delete[] ps_out, ps_in;

	return;
}

// �������W�n����a����f�ʍ��W�n�ւ̕ϊ�
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

// (1) �w�肵���P�ʒ����E�P�ʎ��ԂŋK�i�����ꂽ���Ԕ����l�̎��o���e�X�g�D��v�Z�Ƃ̈�v���m�F�D
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

// (2) �S���������ݒ肳��Ă���ꍇ�̋����̊m�F
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
