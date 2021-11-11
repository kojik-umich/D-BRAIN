#include "pch.h"
#include "B4P_StabIn.h"

#define _FACE_		1e5
#define _CORNEROUT_ 2e5
#define _OPEN_		3e5
#define _CORNERIN_	4e5
#define _EDGEIN_	5e5
#define _EDGEOUT_	6e5

class bal_SnapCageTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	bal_SnapCage CG;
	double BL_r;
	Vector3d BL_Cent;						

	virtual void SetUp() {
		// �����l��50BSWZ�̂��̂�p����D
		B4P_StabIn FI;

		FI.Snap.dout		= 72.2e-3;
		FI.Snap.din			= 62.8e-3;
		FI.Cage.m			= 65.553e-3;

		Vector3d rmg0(3.64545273e-3, 0.0, 0.0);
		for (int i = 0; i < 3; i++)
			FI.Cage.rmg0[i]		= rmg0[i];

		FI.ballnum			= 11;
		FI.Cage.Ix			= 7.5031;
		FI.Cage.Iyz			= 0.022;

		FI.balldia			= 11.1e-3;		// �ʌa[m]
		FI.Snap.R			= 5.74e-3;		// �|�P�b�g���a[m]
		FI.Snap.ropen		= 5.055e-3;		// �J�������a[m]
		FI.Snap.Kface = _FACE_;
		FI.Snap.Kcornerout = _CORNEROUT_;
		FI.Snap.Kopen = _OPEN_;
		FI.Snap.Kcornerin = _CORNERIN_;
		FI.Snap.Kedgein = _EDGEIN_;
		FI.Snap.Kedgeout = _EDGEOUT_;
		FI.Snap.jc			= 0.0;

		this->CG.init(FI);

		// �ێ���􉽒��S�����_�ɗ���悤�ɒ����D
		for (int i = 0; i < 3; i++)
			this->CG.x[i] = -FI.Cage.rmg0[i];
		this->CG.q = Quaterniond::Identity();

		this->BL_r = FI.balldia / 2;
		this->BL_Cent = Vector3d(0.0, 0.0, 67.5e-3 / 2);
	}

	virtual void TearDown() {
	}
};

// (1) �p���W�C��������������
TEST_F(bal_SnapCageTest, init0) {
	// �|�P�b�g0�̍��W
	Vector3d corno0_ex((3.5 - 0.78070335) * 1e-3, -4.61921891e-3, 35.80325148e-3);
	Vector3d corno1_ex((3.5 - 0.78070335) * 1e-3, 4.61921891e-3, 35.80325148e-3);
	Vector3d corni0_ex((3.5 - 0.78070335) * 1e-3, -4.30670957e-3, 31.10325148e-3);
	Vector3d corni1_ex((3.5 - 0.78070335) * 1e-3, 4.30670957e-3, 31.10325148e-3);

	// �|�P�b�g���_�̍��W�͉�]�s���p���Čv�Z���C�v���O�����v�Z�l�Ɣ�r
	VectorXd th = VectorXd(12).setLinSpaced(0, 2 * Numeric::pi);
	for (int i = 0; i < 11; i++) {
		double thi = th[i];
		Matrix3d roti;		// ��]�s��
		roti << 1, 0, 0,
			0, cos(-th[i]), sin(-th[i]),
			0, -sin(-th[i]), cos(-th[i]);
		EXPECT_NEAR((this->CG.PK[i].corno0 - roti * corno0_ex).norm(), 0, 1e-6);
		EXPECT_NEAR((this->CG.PK[i].corno1 - roti * corno1_ex).norm(), 0, 1e-6);
		EXPECT_NEAR((this->CG.PK[i].corni0 - roti * corni0_ex).norm(), 0, 1e-6);
		EXPECT_NEAR((this->CG.PK[i].corni1 - roti * corni1_ex).norm(), 0, 1e-6);
	}


	// �ŏ����������ƊJ���p��(�J���������ŕ]��), ���͂��ׂẴ|�P�b�g�œ������l�ɂȂ��Ă��邱�Ƃ��m�F�D
	for (int i = 0; i < 11; i++){
		EXPECT_NEAR(this->CG.PK[i].zomin, 35.69370222e-3, 1e-6);
		EXPECT_NEAR(this->CG.PK[i].zimin, 30.99370222e-3, 1e-6);
		double h0_ex = 3.5e-3 - 0.78070681e-3;
		EXPECT_NEAR(this->CG.PK[i].h0, h0_ex, 1e-6);

		// cos���͐���������p���Ē�`�ʂ�v�Z
		double cos_gammao_ex = (35.80325148e-3 - 33.75e-3) / 5.055e-3;
		double cos_gammai_ex = (31.10325148e-3 - 33.75e-3) / 5.055e-3;
		
		EXPECT_NEAR(this->CG.PK[i].cos_gammao, cos_gammao_ex, 1e-6);
		EXPECT_NEAR(this->CG.PK[i].cos_gammai, cos_gammai_ex, 1e-6);
	}
	return;
}

// (1) �ʍ��W����́C�o�͂����ڐG�_�̐��ƌ�����]��
// �֐� get_ContactPoint�Ehow_Contact�Ewhere_Contact ���e�X�g�Ώ�
TEST_F(bal_SnapCageTest, get_ContactPoint1) {

	Vector3d BL_x;
	Vector3d return_[_MAX_CONTACT_];
	double k[_MAX_CONTACT_];
	Vector3d dx;
	int nc, ptt[_MAX_CONTACT_];
	Vector3d _pb, e_pb, _pp, e_pp;

	// �P�[�X(A) face 
	dx = Vector3d(0.0, 0.4e-3, 0.0);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	double dot = dx.dot(return_[0] - BL_x) / dx.norm() / (return_[0] - BL_x).norm();
	EXPECT_EQ(nc, 1);										// �ڐG�_�̐���1��
	EXPECT_NEAR(dot, 1, 1e-6);								// �ʂ̃|�P�b�g�������x�N�g���Ɛڋߗʃx�N�g������������
	double La = (return_[0] - this->CG.PK[0].x).norm();
	EXPECT_NEAR(La, this->CG.PK[0].R, 1e-6);				// �|�P�b�g����ɐڐG�_������
	EXPECT_NEAR(k[0], _FACE_, 1e-3);							// �ڐG����

	// �P�[�X(B) �O�փG�b�W������ edgeo
	dx = Vector3d(0.0, 0.4e-3, 0.4e-3);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_EQ(nc, 1);										// �ڐG�_�̐���1��
	double Lb1 = return_[0].norm();
	EXPECT_NEAR(Lb1, this->CG.ro, 1e-6);					// �ێ���O����ɓ_������Ă��邩�m�F
	double Lb2 = (return_[0] - this->CG.PK[0].x).norm();
	EXPECT_NEAR(Lb2, this->CG.PK[0].R, 1e-6);				// �|�P�b�g���ɓ_������Ă��邩�m�F

	// �a����(Z������)����݂��Ƃ��̕��ʃx�N�g�����ʂ̕��ʃx�N�g���Ɠ����ł��邩�m�F
	_pb = return_[0] - this->CG.PK[0].x;
	e_pb = Vector3d(_pb[0], _pb[1], 0).normalized();		// �a��������݂��ڐG�_�̕��ʃx�N�g��(�|�P�b�g���S�)				
	_pp = return_[0] - this->CG.PK[0].x;
	e_pp = Vector3d(_pp[0], _pp[1], 0).normalized();		// �a��������݂��ʂ̕��ʃx�N�g��(�|�P�b�g���S�)
	EXPECT_NEAR((e_pb - e_pp).norm(), 0, 1e-3);
	EXPECT_NEAR(k[0], _EDGEOUT_, 1e-3);							// �ڐG����
		
	// �P�[�X(C) ���փG�b�W������ edgei
	dx = Vector3d(0.0, 0.4e-3, -0.4e-3);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_EQ(nc, 1);										// �ڐG�_�̐���1��
	double Lc1 = return_[0].norm();
	EXPECT_NEAR(Lc1, this->CG.ri, 1e-6);					// �ێ�������̉~����ɓ_������Ă��邩�m�F
	double Lc2 = (return_[0] - this->CG.PK[0].x).norm();
	EXPECT_NEAR(Lc2, this->CG.PK[0].R, 1e-6);				// �|�P�b�g���ɓ_������Ă��邩�m�F

	// �a����(Z������)����݂��Ƃ��̕��ʃx�N�g�����ʂ̕��ʃx�N�g���Ɠ����ł��邩�m�F
	_pb = return_[0] - this->CG.PK[0].x;
	e_pb = Vector3d(_pb[0], _pb[1], 0).normalized();		// �a��������݂��ڐG�_�̕��ʃx�N�g��(�|�P�b�g���S�)				
	_pp = return_[0] - this->CG.PK[0].x;
	e_pp = Vector3d(_pp[0], _pp[1], 0).normalized();		// �a��������݂��ʂ̕��ʃx�N�g��(�|�P�b�g���S�)
	EXPECT_NEAR((e_pb - e_pp).norm(), 0, 1e-3);
	EXPECT_NEAR(k[0], _EDGEIN_, 1e-3);							// �ڐG����

	// �P�[�X(D) �J���������� aperture
	dx = Vector3d(0.4e-3, 0.4e-3, 0.0);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_EQ(nc, 2);										// �ڐG�_�̐���2��
	EXPECT_NEAR(return_[0].x(), this->CG.PK[0].h0, 1);		// ��ʂɓ_������Ă��邩�m�F
	double Ld1 = (return_[0] - this->CG.PK[0].x).norm();
	EXPECT_NEAR(Ld1, this->CG.PK[0].R, 1e-6);				// �|�P�b�g���ɓ_������Ă��邩�m�F

	// ���(X������)����݂��Ƃ��̐ڐG�_�̕��ʃx�N�g�����ʂ̕��ʃx�N�g���Ɠ����ł��邩�m�F�i�����x�N�g���̓|�P�b�g���S��Ɍv�Z�j
	_pb = return_[0] - this->CG.PK[0].x;
	e_pb = Vector3d(0, _pb[1], _pb[2]).normalized();		// ��ʂ���݂��ڐG�_�̕��ʃx�N�g��(�|�P�b�g���S�)				
	_pp = return_[0] - this->CG.PK[0].x;
	e_pp = Vector3d(0, _pp[1], _pp[2]).normalized();		// ��ʂ���݂��ʂ̕��ʃx�N�g��(�|�P�b�g���S�)
	EXPECT_NEAR((e_pb - e_pp).norm(), 0, 1e-3);
	EXPECT_NEAR(k[0], _OPEN_, 1e-3);							// �ڐG����
	
	// �P�[�X(E) �p�O�֑� cornero
	dx = Vector3d(0.4e-3, -0.4e-3, 0.4e-3);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR((return_[0] - this->CG.PK[0].corno0).norm(), 0, 1e-3);
	EXPECT_NEAR((return_[1] - this->CG.PK[0].corno1).norm(), 0, 1e-3);		// �ڐG�_�����_�ƈ�v
	EXPECT_EQ(nc, 2);														// �ڐG�_�̐���2��
	EXPECT_NEAR(k[0], _CORNEROUT_, 1e-3);											// �ڐG����

	// �P�[�X(F) �p���֑� corneri
	dx = Vector3d(0.4e-3, -0.4e-3, -0.4e-3);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR((return_[0] - this->CG.PK[0].corni0).norm(), 0, 1e-3);
	EXPECT_NEAR((return_[1] - this->CG.PK[0].corni1).norm(), 0, 1e-3);		// �ڐG�_�����_�ƈ�v
	EXPECT_EQ(nc, 2);												  		// �ڐG�_�̐���2��
	EXPECT_NEAR(k[0], _CORNERIN_, 1e-3);											// �ڐG����

	// �P�[�X(O) exception�i�����ʒu�j
	dx = Vector3d(0.0, 0.0, 0.0);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_EQ(nc, 0);			// �ڐG�_�Ȃ�

	// �P�[�X(I) exception�i�ێ���􉽒��S�E�|�P�b�g���S�E�{�[�����S�����j
	dx = Vector3d(0.0, 0.0, 1e-3);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_EQ(nc, 0);			// �ڐG�_�Ȃ�

	// �P�[�X(J) exception�i���ƃ|�P�b�g���S-�{�[�����S�����s�j
	dx = Vector3d(1e-3, 0.0, 0.0);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_EQ(nc, 0);			// �ڐG�_�Ȃ�

	return;
};

// (2) ���_�̋��E�l�t�߂̓����蔻��̊m�F�i�ڐG�_�̍��W�̓e�X�g���Ȃ��j
TEST_F(bal_SnapCageTest, get_ContactPoint2) {

	Vector3d BL_x;
	Vector3d return_[_MAX_CONTACT_];
	double k[_MAX_CONTACT_];
	Vector3d dx;
	double dot;
	int nc, ptt[_MAX_CONTACT_];

	// �p�̈ʒu�ɋʈʒu��ݒ�(�O�֑�)
	Vector3d corno0_ex((3.5 - 0.78070335) * 1e-3, -4.61921891e-3, 35.80325148e-3);

	// -x, +z�ł���΃G�b�W������
	dx = Vector3d(-1e-6, 0, 1e-6);
	BL_x = corno0_ex + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _EDGEOUT_, 1e-3);			// �ڐG����
	EXPECT_EQ(nc, 1);						// �ڐG�_�̐���1��

	// +x, +z�ł���Ίp������
	dx = Vector3d(1e-6, 0, 1e-6);
	BL_x = corno0_ex + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _CORNEROUT_, 1e-3);			// �ڐG����
	EXPECT_EQ(nc, 2);						// �ڐG�_�̐���2��

	// +x, -z�ł���ΊJ����������
	dx = Vector3d(1e-6, 0, -1e-6);
	BL_x = corno0_ex + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _OPEN_, 1e-3);			// �ڐG����
	EXPECT_EQ(nc, 2);						// �ڐG�_�̐���2��


	// �p�̈ʒu�ɋʈʒu��ݒ�(���֑�)
	Vector3d corni0_ex((3.5 - 0.78070335) * 1e-3, -4.30670957e-3, 31.10325148e-3);

	// -x, -z�ł���΃G�b�W������
	dx = Vector3d(-1e-6, 0, -1e-6);
	BL_x = corni0_ex + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _EDGEIN_, 1e-3);			// �ڐG����
	EXPECT_EQ(nc, 1);						// �ڐG�_�̐���1��

	// +x, -z�ł���Ίp������
	dx = Vector3d(1e-6, 0, -1e-6);
	BL_x = corni0_ex + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _CORNERIN_, 1e-3);			// �ڐG����
	EXPECT_EQ(nc, 2);						// �ڐG�_�̐���2��

	// +x, +z�ł���ΊJ����������
	dx = Vector3d(1e-6, 0, 1e-6);
	BL_x = corni0_ex + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _OPEN_, 1e-3);			// �ڐG����
	EXPECT_EQ(nc, 2);						// �ڐG�_�̐���2��

	return;
}

// (3) �J�����̋��E�l�t�߂̓����蔻��̊m�F�i�ڐG�_�̍��W�̓e�X�g���Ȃ��j
TEST_F(bal_SnapCageTest, get_ContactPoint3) {
	
	Vector3d BL_x;
	Vector3d return_[_MAX_CONTACT_];
	double k[_MAX_CONTACT_];
	Vector3d dx;
	double dot;
	int nc, ptt[_MAX_CONTACT_];

	// �J�����̍��W
	Vector3d _x = Vector3d(this->CG.PK[0].h0, this->CG.PK[0].ropen, 0) 
		+ this->CG.PK[0].x;

	// -x�����ɓ����Ɩʓ�����
	dx = Vector3d(-1e-6, 0, 0);
	BL_x = _x + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _FACE_, 1e-3);			// �ڐG����
	EXPECT_EQ(nc, 1);						// �ڐG�_�̐���1��

	// +x�����ɓ����ƊJ����������
	dx = Vector3d(1e-6, 0, 0);
	BL_x = _x + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _OPEN_, 1e-3);			// �ڐG����
	EXPECT_EQ(nc, 2);						// �ڐG�_�̐���2��
	
	return;
}

// (4) ���E�l���ׂ��t�߁i�ʁ��G�b�W���R�[�i�[���J�������ʁj�̐ڐG�_�̈ړ������C���̋ߖT�ŘA���ɂȂ邩�̊m�F�D
TEST_F(bal_SnapCageTest, get_ContactPoint4) {

	double r = 1e-4;
	int n = 50;
	VectorXd th = VectorXd(n).setLinSpaced(0, 2* Numeric::pi);
	double dth = abs(th[1] - th[0]);
	double dx = r * dth;	// ���悻 1.282e-05 m �� 12 um

	// �܂��́C���R�[�i�[����ŉ~�`�ɋʈʒu�𓮂����C�ڐG�_�̈ړ����m�F����D�D
	Vector3d corni0 = this->CG.PK[0].corni0;
	Vector3d ci0 = this->CG.PK[0].corni0 - this->CG.PK[0].x;
	Vector3d axi = Vector3d::UnitX().cross(ci0).normalized();
	Vector3d ayi = axi.cross(ci0).normalized();
	int ptt[_MAX_CONTACT_];
	Vector3d xc0;
	for (size_t i = 0; i < n; i++) {
		Vector3d BL_x = corni0 + r * cos(th[i]) * axi + r * sin(th[i]) * ayi;

		Vector3d xc[_MAX_CONTACT_];		double k[_MAX_CONTACT_];
		int nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, xc, k, ptt);
		if (i > 0)
			EXPECT_LT((xc0 - xc[0]).norm(), 3 * dx);	// ������1*dx�ɕς���Ƃ������e�X�g��ʂ�Ȃ����̂��o�Ă��܂��D���̃R�����g�ƃZ�b�g�ł��Ƃ�蕪����₷���ł��D�����āC������ƍl���Ă݂Ă݂Ă���������
		xc0 = xc[0];
		// cout << k[0] << endl;	// �������o�͂����邱�ƂŁC�ڐG���肪�ǂ��Ȃ��Ă���̂��m�F�ł��܂��D����m�F���ɂ͂��̃R�����g���O���Ď��s���Ă��������D
	}

	// ���l�ɁC�O�R�[�i�[����D
	Vector3d corno0 = this->CG.PK[0].corno0;
	Vector3d co0 = this->CG.PK[0].corno0 - this->CG.PK[0].x;
	Vector3d axo = Vector3d::UnitX().cross(co0).normalized();
	Vector3d ayo = axo.cross(co0).normalized();

	for (size_t i = 0; i < n; i++) {
		Vector3d BL_x = corno0 + r * cos(th[i]) * axo + r * sin(th[i]) * ayo;

		Vector3d xc[_MAX_CONTACT_];		double k[_MAX_CONTACT_];
		int nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, xc, k, ptt);
		if (i > 0)
			EXPECT_LT((xc0 - xc[0]).norm(), 3 * dx);
		xc0 = xc[0];
		// cout << k[0] << endl;	// �������o�͂����邱�ƂŁC�ڐG���肪�ǂ��Ȃ��Ă���̂��m�F�ł��܂��D����m�F���ɂ͂��̃R�����g���O���Ď��s���Ă��������D
	}
	return;
}

