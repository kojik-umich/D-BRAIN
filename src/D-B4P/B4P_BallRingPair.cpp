#include "B4P_BallRingPair.h"

// �ʁE���ցi�O�ցj�I�u�W�F�N�g�����L������
void B4P_BallRingPair::link(Ball*BL, B4P_Ring*RG) {
	this->BL = BL;
	this->RG = RG;

	return;
}

// ���݂̏�Ԃł̋ʂ̐ڐG���肨��сC���̂Ƃ��̐ڐG��Ԃ��v�Z
// �Q�l "01. �~�ʂƋ��̂̒e���ڐG�A���S���Y��.pptx"
bool B4P_BallRingPair::how_Contact
(							// out:	[bit]	�ڐG���Ă����true�C���Ă��Ȃ����false��Ԃ��D
	int i,					// in:	[-]		�ڐG�a�ԍ��D0 or 1�D
	const Vector3d&bl_x,	// in:	[m]		�{�[���ʒu�D�i�����O���W�n�j
	Vector3d&er,			// out:	[-]:	�����O���S���猩���{�[���ʑ������x�N�g���Dx����=0�D(�����O���W�n�)
	Vector3d&eg,			// out:	[-]:	�a���S���猩���{�[�������x�N�g���D�{�[���̎󂯂�ڐG�׏d�����̋t�����D�i�����O���W�n��j
	double&dx				// out:	[-]		�e���ڋߗʁD���ŐڐG�D
) {

	// �܂�yz���ʂŌ����{�[���̕����x�N�g�� er (�����O���W�n�)�𓱏o�D
	double r = sqrt(Numeric::Square(bl_x[1]) + Numeric::Square(bl_x[2]));
	er = Vector3d(0.0, bl_x[1] / r, bl_x[2] / r);

	// ���ڂ���a���S�_�i�����O���W�n��j������D
	Vector3d x_trs = this->RG->GV[i].Rr *er;
	x_trs[0] = this->RG->GV[i].Rx;

	// �a���S�_����{�[�����S�Ɍ������x�N�g���̎Z�o�D
	Vector3d gb = bl_x - x_trs;
	double   gb_norm = gb.norm();

	// �a���S�_�ƃ{�[�����S����v���Ă�����ڂ��Ă��Ȃ����Ƃɂ���D�i0���Z�h�~�j
	if (gb_norm == 0)
		return false;

	eg = gb / gb_norm;
	// �􉽊w����{�[���H�����ݗʂ̎Z�o�D
	dx = gb_norm - this->RG->GV[i].r + this->BL->r;

	// �{�[�����S���g�[���X����o�Ă�����ڂ��Ă��Ȃ����Ƃɂ���D
	if (gb_norm > this->RG->GV[i].r) {
		cout << "�ʂ����O�ւ̋O������O��Ă��܂��D�����Ɍv�Z���~���Ă��������D" << endl;
		return false;
	}

	// �ڋߗʂ����l�̏ꍇ�ڐG���Ă��Ȃ��Ddx=0�̂Ƃ����ȍ~�̌v�Z��0�����������邽�߁C�ڐG�Ȃ��Ƃ���D
	if (dx <= 0)
		return false;

	return true;
}


// Hertz�̌v�Z����ڐG������Ԃ����\�b�h�i�����Ȃ��j
// �Q�l "01. �~�ʂƋ��̂̒e���ڐG�A���S���Y��.pptx" STEP 4 ~ 7
void B4P_BallRingPair::calc_Hertz
(							// out:	[-]		�Ȃ�
	int i,					// in:	[-]		�ڐG�a�ԍ��D0 or 1�D
	const Vector3d&bl_x,	// in:	[m]		�{�[���ʒu�D�i�����O���W�n�j
	const Vector3d&er,		// in:	[-]:	�����O���S���猩���{�[���ʑ������x�N�g���Dx����=0�D(�����O���W�n)
	const Vector3d&eg,		// in:	[-]:	�a���S���猩���{�[�������x�N�g���D�{�[���̎󂯂�ڐG�׏d�����̋t�����D(�����O���W�n)
	double dx,				// in:	[-]		�e���ڋߗʁD���ŐڐG�D
	double&Rx,				// out:	[m]:	�]��������ȗ����a�D
	double&Ry,				// out:	[m]:	�������ȗ����a�D
	Vector3d&p,				// out:	[m]:	�ڐG�_�ʒu�i�������W�n�j�D
	double&cos_alp,			// out:	[-]:	�ڐG�p��cos�D
	double&sin_alp,			// out:	[-]:	�ڐG�p��sin�D
	double&a,				// out:	[m]:	�ڐG�ȉ~���a�D
	double&b,				// out:	[m]:	�ڐG�ȉ~�Z�a�D
	double&k				// out:	[N/x^1.5]:	����`�����D
) {
	// �ڐG�p alp �� cos �̎Z�o�D�i�i�b�g���W�A��������=0�x�C�V���t�g���W�A��������=180�x�j
	Vector3d ea(1.0, 0.0, 0.0);
	cos_alp = er.dot(eg);
	sin_alp = ea.dot(eg);

	// �ȗ��̎Z�o�D�iBrewHamrock �̋ߎ����Ŏg�p�j
	double rho0 = this->BL->r_inv;			// �{�[���̎������ȗ��D
	double rho1 = this->BL->r_inv;			// �{�[���̌a�����ȗ��D
	double rho2 = -cos_alp / (this->RG->GV[i].Rr + this->RG->GV[i].r * cos_alp);	// �ւ̓]��������ȗ��D
	double rho3 = -this->RG->GV[i].r_inv;	// �ւ̎������ȗ��D

	Rx = 1. / (rho0 + rho2);
	Ry = 1. / (rho1 + rho3);

	// BrewHamrock �̋ߎ������獄���C�ڐG�ȉ~�����߂�D 
	this->HZ->calc(Rx, Ry, dx, this->E, k, a, b);

	// ��r�p�D�ȑO�܂ł̎d�l�ł��D�m�F���I�������ȉ�2�s�͏����Ă��������D
	//double Rm = 2 * this->BL->r * this->RG->GV[i].r / (this->BL->r + this->RG->GV[i].r);	// �ڐG�ʂ̋ȗ����a
	//Vector3d p__ = bl_x + (this->BL->r + a * a * 0.5 * (1 / Rm - 1 / this->BL->r)) * eg;

	// �ڐG�_�ʒup�̎Z�o�D�i�����O���W�n�ˊ������W�n�D�j
	double x0 = Numeric::EffectiveCenter(this->BL->r, this->RG->GV[i].r, dx);
	Vector3d p_ = bl_x + x0 * eg;
	p = this->RG->to_inecoord(p_);

	return;
}

// Hertz�̌v�Z����C�_���p�ɂ�錸�����l�������׏d��Ԃ����\�b�h�D
double B4P_BallRingPair::calc_DynamicHertz
(							// out:	[N]		�ڐG�׏d�i�X�J���[�j
	int i,					// in:	[-]		�ڐG�a�ԍ��D0 or 1�D
	const Vector3d&bl_x,	// in:	[m]		�{�[���ʒu�D�i�����O���W�n�j
	const Vector3d&bl_v,	// in:	[m/s]	�{�[�����x�D�i�����O���W�n�j
	const Vector3d&er,		// in:	[-]:	�����O���S���猩���{�[���ʑ������x�N�g���Dx����=0�D(�����O���W�n)
	const Vector3d&eg,		// in:	[-]:	�a���S���猩���{�[�������x�N�g���D�{�[���̎󂯂�ڐG�׏d�����̋t�����D(�����O���W�n)
	double dx,				// in:	[-]		�e���ڋߗʁD���ŐڐG�D
	double&Rx,				// out:	[m]:	�]��������ȗ����a�D
	double&Ry,				// out:	[m]:	�������ȗ����a�D
	Vector3d&p,				// out:	[m]:	�ڐG�_�ʒu�i�������W�n�j�D
	double&cos_alp,			// out:	[-]:	�ڐG�p��cos�D
	double&sin_alp,			// out:	[-]:	�ڐG�p��sin�D
	double&a,				// out:	[m]:	�ڐG�ȉ~���a�D
	double&b				// out:	[m]:	�ڐG�ȉ~�Z�a�D
) {
	// Hertz�ڐG�������v�Z
	double k;
	this->calc_Hertz(i, bl_x, er, eg, dx, Rx, Ry, p, cos_alp, sin_alp, a, b, k);

	// �ʑ��x�C�ڐG�ʐ��������i�����O���W�n�ƃ����O���W�n�̓��ρj
	double bl_vn = bl_v.dot(eg);

	// ����`�̌�����(�X�J���[)
	double F_norm_ = this->DF->calc(k, this->zeta, this->m, bl_vn, dx);
	return F_norm_;
}

// ���薀�C�����߂郁�\�b�h�C�g���N�͋O���ւɂ�������̂Ƌʂɂ�������̂𗼕��v�Z
void B4P_BallRingPair::calc_Sliding
(
	int i,				// in:	[-]		�ڐG�a�ԍ��D0 or 1�D
	const Vector3d&p,	// in:	[m]:	�ڐG�_�ʒu�i�������W�n�j�D
	double a,			// in:	[m]:	�ڐG�ȉ~���a�D
	double Pmax,		// in:	[Pa]:	�ő�ʈ��D
	double F_norm,		// in:	[N]:	�w���c�ڐG�׏d�D
	double fratio,		// in:	[-]:	�����ڐG�����D
	Vector3d&Fbs,		// out:	[N]:	�{�[���ɂ����銊�薀�C�́i�������W�n�j
	Vector3d&Tbs,		// out:	[Nm]:	�{�[���ɂ����銊�薀�C�g���N�i�������W�n�j
	Vector3d&Tis		// out:	[Nm]:	����(�O��)�ɂ����銊�薀�C�g���N�i�������W�n�j
) {
	// �ڐG�ȉ~���X���C�X���Ă��ꂼ��̊��薀�C�͂�ώZ����D
	Fbs = Tbs = Tis = Vector3d::Zero();
	this->RG->calc_slice(p, a, this->GV[i].msmax, i, this->BL->r, this->GV[i].ps);

	for (int j = 0; j < this->GV[i].msmax; j++) {
		Vector3d us_, ur_;
		double f_arr = F_norm * this->GV[i].ratio_slice[j];
		this->get_us_ur(this->GV[i].ps[j], us_, ur_);
		double us_norm = us_.norm(), ur_norm = ur_.norm();
		double mu_tr = this->TR->calc(this->lb_eta, Pmax, us_norm, ur_norm);
		double mu_cl = this->CL->calc(this->mu, us_norm, ur_norm, this->clmb_s);
		double Fs_norm = mu_tr * fratio * f_arr + mu_cl * (1.0 - fratio) * f_arr;
		Vector3d Fs_ = Fs_norm * -us_ / us_norm;
		Vector3d Tbs_ = this->BL->calc_Torque(this->GV[i].ps[j], Fs_);
		Vector3d Tis_ = this->RG->calc_Torque(this->GV[i].ps[j], -Fs_);
		Fbs += Fs_;
		Tbs += Tbs_;
		Tis += Tis_;
		this->save_Slice(i, j, f_arr, Fs_, Tbs_, mu_cl, mu_tr, us_, this->GV[i].ps[j]);
	}
	return;
}

// �v�Z���ʂ������o�ϐ��Ɋi�[���郁�\�b�h�D
void B4P_BallRingPair::save(bool c, int i, double cos_alp, double sin_alp, double dx, double a, double b, const Vector3d&p, const Vector3d&Fn, double Pmax, const Vector3d&us, const Vector3d&ur, const Vector3d&Tr, double lambda, double fratio, const Vector3d&Fs, const Vector3d&Ts, const Vector3d&Fb, const Vector3d&Tb, const Vector3d&Ti
) {
	this->GV[i].cos_alp = cos_alp;
	this->GV[i].sin_alp = sin_alp;
	this->GV[i].dx = dx;
	this->GV[i].a = a;
	this->GV[i].b = b;
	this->GV[i].p = p;
	this->GV[i].Fn = Fn;
	this->GV[i].Pmax = Pmax;
	this->GV[i].us = us;
	this->GV[i].ur = ur;
	this->GV[i].Tr = Tr;
	this->GV[i].lambda = lambda;
	this->GV[i].fratio = fratio;
	this->GV[i].Fs = Fs;
	this->GV[i].Ts = Ts;
	this->GV[i].Fb = Fb;
	this->GV[i].Tb = Tb;
	this->GV[i].Ti = Ti;

	return;
}
// �X���C�X�Ђ̌��ʂ�0���i�[
void B4P_BallRingPair::init_Sliceparam(int i) {
	for (int j = 0; j < this->GV[i].msmax; j++) {
		this->GV[i].SL[j].f_arr = 0;
		this->GV[i].SL[j].fs = Vector3d::Zero();	// ���薀�C��[N](�������W�n)
		this->GV[i].SL[j].ts = Vector3d::Zero();	// ���薀�C�g���N[Nm](�������W�n)
		this->GV[i].SL[j].mu_cl = 0;				// �N�[�������C�W��[-]
		this->GV[i].SL[j].mu_tr = 0;				// �g���N�V�����W��[-]
		this->GV[i].SL[j].us = Vector3d::Zero();	// ���葬�x[m/s](�������W�n)
		this->GV[i].SL[j].ps = Vector3d::Zero();	// ���葬�x[m/s](�������W�n)			
	}
}



// �X���C�X�v�Z���ʂ������o�ϐ��Ɋi�[
void B4P_BallRingPair::save_Slice(
	int i,			// in: �a�ԍ�
	int j,			// in: �X���C�X�ԍ�
	double f_arr,	// in: �e�X���C�X�̐ڐG�׏d[N]
	const Vector3d& Fs_,	// in: ���薀�C�׏d[N]
	const Vector3d& Ts_,	// in: ���薀�C�g���N[Nm]
	double mu_cl,	// in: �N�[�������C�W��[-]
	double mu_tr,	// in: �g���N�V�����W��[-]
	const Vector3d& us_,	// in: ���葬�x[m/s]
	const Vector3d& ps_		// in: �X���C�X�В��S[m]
) {
	this->GV[i].SL[j].f_arr = f_arr;
	this->GV[i].SL[j].fs = Fs_;			// ���薀�C��[N](�������W�n)
	this->GV[i].SL[j].ts = Ts_;			// ���薀�C�g���N[Nm](�������W�n)
	this->GV[i].SL[j].mu_cl = mu_cl;		// �N�[�������C�W��[-]
	this->GV[i].SL[j].mu_tr = mu_tr;		// �g���N�V�����W��[-]
	this->GV[i].SL[j].us = us_;			// ���葬�x[m/s](�������W�n)
	this->GV[i].SL[j].ps = ps_;		// �X���C�X�Ђ̑��Έʒu[m](�������W�n)
	return;
}

// ���݂̃{�[��-�����O�̕ψʂ���C���݂��ɂ�����׏d���Z�o����D
void B4P_BallRingPair::calc_force(
	Vector3d&Fbi,	// out:	[N]:	�{�[���ɂ�����S�Ắ̗i�������W�n�j
	Vector3d&Tbi,	// out:	[Nm]:	�{�[���ɂ�����S�Ẵg���N�i�������W�n�j
	Vector3d&Fib, 	// out:	[N]:	�����O�ɂ�����S�Ắ̗i�������W�n�j
	Vector3d&Tib	// out:	[Nm]:	�����O�ɂ�����S�Ẵg���N�i�������W�n�j
) {
	Vector3d bl_x = this->RG->to_mycoord(this->BL->x);
	Vector3d bl_v = this->RG->to_myvelocity(this->BL->v);
	Fbi = Tbi = Fib = Tib = Vector3d::Zero();

	// groove �� ball �̐ڐG�v�Z�D
	for (int i = 0; i < 2; i++) {

		// �K�v�Ȗ߂�l�̊m�ہD
		double dx, Rx, Ry, cos_alp, sin_alp, a, b, Pmax, lambda, fratio, k;
		dx = Rx = Ry = sin_alp = a = b = Pmax = lambda = fratio = k = 0;
		cos_alp = 1.0;
		Vector3d er, eg, p, Fn, ur, us, ur3, us3, Tbr, Fbs, Tbs, Tis, Fb, Tb, Ti;
		er = eg = p = Fn = ur = us = Tbr = Fbs = Tbs = Tis = Fb = Tb = Ti = Vector3d::Zero();

		// �ڂ��Ă��Ă��邩����D
		bool c = this->how_Contact(i, bl_x, er, eg, dx);
		if (c) {
			// �ڐG�׏d�E�]���̖ʈ��𓱏o�i�������W�n�j
			double F_norm = this->calc_DynamicHertz(i, bl_x, bl_v, er, eg, dx, Rx, Ry, p, cos_alp, sin_alp, a, b);
			Fn = F_norm * this->RG->to_inevector(-eg);
			Pmax = 1.5 * F_norm / (Numeric::pi * a * b);

			//// ���ڂ���a���S�_�i�����O���W�n��j������D
			//Vector3d x_trs = this->RG->GV[i].Rr *er;
			//x_trs[0] = this->RG->GV[i].Rx;
			//// �a���S�_����{�[�����S�Ɍ������x�N�g���̎Z�o�D
			//Vector3d gb = bl_x - x_trs;
			//
			//Vector3d _bl_x = Vector3d(0, bl_x[1], bl_x[2]);
			//Vector3d _bl_v = Vector3d(0, bl_v[1], bl_v[2]);
			//Vector3d _dpdt = bl_v + this->BL->r * (_bl_v - this->RG->GV[i].Rr * _bl_v / _bl_x.norm()) / gb.norm();
			//Vector3d dpdt = this->RG->to_inevelocity(_dpdt);
			// �]���葬�x�C���葬�x�̎Z�o�D�i�������W�n�j
			//this->get_us_ur2(p, dpdt, us, ur);
			//this->get_us_ur3(p, er, us, ur);
			//Vector3d ey = Vector3d(0, 1, 0);	// �a�������W�n��Y�����P�ʃx�N�g��
			//Vector3d thXZ = this->RG->ine_to_XZcoord(this->BL->x);
			//Matrix3d R_th = this->RG->get_xyz2XYZ(thXZ[0]);
			//Matrix3d R_th_ = R_th.inverse();
			//Vector3d ey_ = R_th_ * ey;
			//Vector3d ey__ = this->RG->to_inevector(ey_);
			//double dot = ey_.dot(er);
			//Vector3d ex = Vector3d(1, 0, 0);
			//Vector3d e_th = this->RG->to_inevector(ex.cross(er));
			//Vector3d e_th_ = e_th / e_th.norm();
			//Vector3d ur_new = ey__ * (ey__.dot(ur));
			//double dot2 = ey__.dot(this->BL->x - this->RG->x);
			
			this->get_us_ur(p, us, ur);
			double ur_norm = ur.norm();

			// �{�[���ɂ����銊�薀�C�R�͂̎Z�o�D
			double h = this->FT->calc(F_norm, this->E, Rx, Ry, this->lb_alpha, this->lb_eta, ur_norm, this->lm, b);
			h *= Tribology::ErtelGrubin(this->lb_eta, this->lb_beta, this->lb_k, ur_norm);
			lambda = h / this->sigma;
			fratio = Tribology::ForceRatio(lambda);
			this->calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);


			// �{�[���ɂ�����]���薀�C�R�͂̎Z�o�D�i�������W�n�j
			double Trf = this->RR->calc(Rx, Ry, this->BL->D, a, b, F_norm, this->E, ur_norm, fratio, this->lb_eta, this->lb_alpha, this->lb_beta, this->lb_k, this->lm);
			double Trh = this->HY->calc(this->fh, b, F_norm);
			Tbr = (Trf + Trh) * this->BL->calc_TorqueDirection(p, -ur);

			Fb = Fn + Fbs;
			Tb = Tbr + Tbs;
			Ti = this->RG->calc_Torque(p, -Fn) + Tis - Tbr;
		}
		else {
			// �X���C�X�Ђ̌��ʂ����ׂ�0�ɂ���
			this->init_Sliceparam(i);
		}

		// �v�Z���ʂ������o�ϐ��ɕۑ�����D
		this->save(c, i, cos_alp, sin_alp, dx, a, b, p, Fn, Pmax, us, ur, Tbr, lambda, fratio, Fbs, Tbs, Fb, Tb, Ti);

		// �a0�ƍa1�̍��v���o�͂���D
		Fbi += Fb;		// Fbi : �{�[��(b)�������O(i)����󂯂�́D
		Tbi += Tb;		// Tbi : �{�[��(b)�������O(i)����󂯂郂�[�����g�D
		Fib += -Fb;		// Fib : �����O(i)���{�[��(b)����󂯂�́D
		Tib += Ti;		// Tib : �����O(i)���{�[��(b)����󂯂郂�[�����g�D
	}
	return;
}

// �v�Z���ʂ������o�ϐ�����o�͗p�ϐ��ɓn��
void B4P_BallRingPair::write(
	B4P_Out::BallRingPair&BRP	// out: (�\����): ��-���O�֐ڐG�v�Z�o��
) {

	this->thXZ = this->RG->ine_to_XZcoord(this->BL->x);
	BRP.th = this->thXZ[0];
	BRP.X = this->thXZ[1];
	BRP.Z = this->thXZ[2];
	Matrix3d xyz2XYZ = this->RG->get_xyz2XYZ(BRP.th);

	for (int ig = 0; ig < 2; ig++) {
		Vector3d Fr = this->GV[ig].Tr.norm() * this->BL->r_inv * -this->GV[ig].ur.normalized();
		for (int j = 0; j < 3; j++) {
			BRP.GV[ig].p[j] = this->GV[ig].p[j];
			BRP.GV[ig].Fn[j] = this->GV[ig].Fn[j];
			BRP.GV[ig].Fs[j] = this->GV[ig].Fs[j];
			BRP.GV[ig].Ts[j] = this->GV[ig].Ts[j];
			BRP.GV[ig].Fr[j] = Fr[j];
			BRP.GV[ig].Tr[j] = this->GV[ig].Tr[j];
			BRP.GV[ig].us[j] = this->GV[ig].us[j];
			BRP.GV[ig].ur[j] = this->GV[ig].ur[j];
		}
		Vector3d Fn_ = this->RG->ine_to_XZvector(this->GV[ig].Fn, xyz2XYZ);
		Vector3d Fs_ = this->RG->ine_to_XZvector(this->GV[ig].Fs, xyz2XYZ);
		Vector3d Ts_ = this->RG->ine_to_XZvector(this->GV[ig].Ts, xyz2XYZ);
		Vector3d Fr_ = this->RG->ine_to_XZvector(Fr, xyz2XYZ);
		Vector3d Tr_ = this->RG->ine_to_XZvector(this->GV[ig].Tr, xyz2XYZ);
		Vector3d us_ = this->RG->ine_to_XZvector(this->GV[ig].us, xyz2XYZ);
		Vector3d ur_ = this->RG->ine_to_XZvector(this->GV[ig].ur, xyz2XYZ);
		for (int j = 0; j < 3; j++) {
			BRP.GV[ig].Fn_[j] = Fn_[j];
			BRP.GV[ig].Fs_[j] = Fs_[j];
			BRP.GV[ig].Ts_[j] = Ts_[j];
			BRP.GV[ig].Fr_[j] = Fr_[j];
			BRP.GV[ig].Tr_[j] = Tr_[j];
			BRP.GV[ig].us_[j] = us_[j];
			BRP.GV[ig].ur_[j] = ur_[j];
		}
		Vector3d p_ = this->RG->ine_to_XZcoord(this->GV[ig].p);
		for (int j = 0; j < 2; j++)
			BRP.GV[ig].p_[j] = p_[j + 1];

		BRP.GV[ig].lambda = this->GV[ig].lambda;
		BRP.GV[ig].fratio = this->GV[ig].fratio;
		BRP.GV[ig].phi = this->ContactAngle(this->GV[ig].cos_alp, this->GV[ig].sin_alp);
		BRP.GV[ig].a = this->GV[ig].a;
		BRP.GV[ig].b = this->GV[ig].b;
		BRP.GV[ig].dx = this->GV[ig].dx;
		BRP.GV[ig].Pmax = this->GV[ig].Pmax;

		this->write_slice(ig, xyz2XYZ, BRP);
	}

	return;
}

// �X���C�X�Ђ̌��ʂ��o�͗p�\���̂Ɋi�[
void B4P_BallRingPair::write_slice(
	int ig,						// in:	[-]: �a�ԍ�
	Matrix3d xyz2XYZ,			// in:	[-]: �a�������W�n�ϊ��s��
	B4P_Out::BallRingPair&BRP	// out:    : ��-���O�֐ڐG�v�Z�o��
) {
	for (int j = 0; j < this->GV[ig].msmax; j++) {
		BRP.GV[ig].SL[j].f_arr = this->GV[ig].SL[j].f_arr;
		BRP.GV[ig].SL[j].mu_cl = this->GV[ig].SL[j].mu_cl;
		BRP.GV[ig].SL[j].mu_tr = this->GV[ig].SL[j].mu_tr;
		Vector3d fs_ = this->RG->ine_to_XZvector(this->GV[ig].SL[j].fs, xyz2XYZ);
		Vector3d ts_ = this->RG->ine_to_XZvector(this->GV[ig].SL[j].ts, xyz2XYZ);
		Vector3d us_ = this->RG->ine_to_XZvector(this->GV[ig].SL[j].us, xyz2XYZ);
		Vector3d p_ = this->RG->ine_to_XZvector(this->GV[ig].SL[j].ps, xyz2XYZ);

		for (int k = 0; k < 3; k++) {
			BRP.GV[ig].SL[j].fs_[k] = fs_[k];
			BRP.GV[ig].SL[j].ts_[k] = ts_[k];
			BRP.GV[ig].SL[j].us_[k] = us_[k];
			BRP.GV[ig].SL[j].ps_[k] = p_[k];
		}
	}
	return;
}
double B4P_BallOuterRingPair::ContactAngle(double cos_alp, double sin_alp) {
	double alp = atan2(sin_alp, cos_alp);
	return alp;
}

double B4P_BallInnerRingPair::ContactAngle(double cos_alp, double sin_alp) {
	double alp = atan2(sin_alp, -cos_alp);
	return alp;
}

// �]���葬�x���Z�o���郁�\�b�h�i�������W�n�j
Vector3d B4P_BallRingPair::get_ur(
	const Vector3d&p		// �ڐG�_�ʒu�D�i�������W�n�j
) {
	Vector3d rwb = this->BL->w.cross(p - this->BL->x);
	Vector3d rwr = this->RG->w.cross(p - this->RG->x);
	Vector3d ur = 0.5 * (3 * this->RG->v - this->BL->v + rwb + rwr);
	Vector3d ur_ = ur;// this->BL->remove_Normal(p, ur);
	return ur_;
}

// ���葬�x���Z�o���郁�\�b�h�i�������W�n�j
Vector3d B4P_BallRingPair::get_us(
	const Vector3d&p		// �ڐG�_�ʒu�D�i�������W�n�j
) {
	Vector3d ub = this->BL->surface_velocity(p);
	Vector3d ui = this->RG->surface_velocity(p);
	Vector3d us = ub - ui;
	Vector3d us_ = us;// this->BL->remove_Normal(p, us);
	return us_;
}

// ���葬�x�Ɠ]���葬�x���Z�o���郁�\�b�h�i�������W�n�j�D���ʂ̕ϐ����������ߍ������̂��߂Ɉ�̂Ƃ����D
void B4P_BallRingPair::get_us_ur(
	const Vector3d&p,		// �ڐG�_�ʒu�D�i�������W�n�j
	Vector3d & us,			// out:	[-]:	���葬�x�D�i�������W�n�j
	Vector3d & ur			// out:	[-]:	�]���葬�x�D�i�������W�n�j
) {
	Vector3d dx = p - this->BL->x;
	Vector3d e = dx.normalized();
	Vector3d rwb = this->BL->w.cross(dx);
	Vector3d rwr = this->RG->w.cross(p - this->RG->x);
	Vector3d ur_ = 0.5 * (3 * this->RG->v - this->BL->v + rwb + rwr);
	ur = ur_;// -e//* e.dot(ur_);

	Vector3d ub = this->BL->v + rwb;
	Vector3d ui = this->RG->v + rwr;
	Vector3d us_ = ub - ui;
	us = us_;// -e; //* e.dot(us_);

	return;
}



// �]���葬�x���Z�o���郁�\�b�h�i�������W�n�j(��2�āC�ڐG�_���x����萳�m�Ɍv�Z���Ă��邪�C���Z�ʂ�����)
void B4P_BallRingPair::get_us_ur2(
	const Vector3d&p,		// in:	[-]:	�ڐG�_�ʒu�D�i�������W�n�j
	const Vector3d&er,		// in:	[-]:	�����O���S���猩���{�[���ʑ������x�N�g���Dx����=0�D
	Vector3d & us,			// out:	[-]:	���葬�x�D�i�������W�n�j
	Vector3d & ur			// out:	[-]:	�]���葬�x�D�i�������W�n�j
) {
	Vector3d dx = p - this->BL->x;
	Vector3d e = dx.normalized();
	Vector3d rwb = this->BL->w.cross(dx);
	Vector3d rwr = this->RG->w.cross(p - this->RG->x);
	Vector3d bl_v = this->BL->v - this->RG->v;
	Vector3d e_a = this->RG->get_ax();
	Vector3d e_th = er.cross(e_a);				// �������x�N�g���i�a������12���̌����̂Ƃ�3���̌����j
	Vector3d bl_vth = bl_v.dot(e_th)*e_th;

	Vector3d RB = this->BL->x - this->RG->x;
	// (���������x�ɂ����������Ȃ���΂Ȃ�Ȃ����C�܂��l���ł��Ă��Ȃ�)
	Vector3d p_v = this->BL->v - this->RG->v + (bl_vth * dx.dot(er) / RB.dot(er));


	Vector3d ur_ = 0.5 * (this->BL->v + this->RG->v - 2 * p_v + rwb + rwr);
	ur = ur_;// -e * e.dot(ur_);

	Vector3d ub = this->BL->v + rwb;
	Vector3d ui = this->RG->v + rwr;
	Vector3d us_ = ub - ui;
	us = us_;// -e * e.dot(us_);

	return;
}

// B4P_In ����ǂݍ���ŏ����ݒ�D
void B4P_BallOuterRingPair::init(const B4P_In&FI) {

	this->mu = FI.BOP.mu;		// ���C�W��
	this->zeta = FI.BOP.dzeta;		// ������

	this->init_Lubrication(FI.LB);
	this->init_Slice(FI.msmax);
	this->init_Tribology(FI.TB);

	return;
}

// B4P_In ����ǂݍ���ŏ����ݒ�D
void B4P_BallInnerRingPair::init(const B4P_In&FI) {

	this->mu = FI.BIP.mu;		//���C�W��
	this->zeta = FI.BIP.dzeta;		//������

	this->init_Lubrication(FI.LB);
	this->init_Slice(FI.msmax);
	this->init_Tribology(FI.TB);

	return;
}

// ���͂ɂ���ė��_�����Z�b�g����D
void B4P_BallRingPair::init_Tribology(const B4P_In::Tribology&FI) {

	this->m = Tribology::ReducedMass(this->BL->m, this->RG->m);									// ���Z����
	this->E = Tribology::ReducedYoung(this->BL->E, this->BL->nu, this->RG->E, this->RG->nu);	// ���������O��
	this->sigma = Tribology::CompositeRoughness(this->BL->sigmap, this->RG->sigmap);			// �����e��

	switch (FI.rollingresistance) {
	case B4P_In::Tribology::RollingResistance::RollingResistanceNothing:
		this->RR = new Tribology::RollingResistanceNothing();
		break;
	case B4P_In::Tribology::RollingResistance::Aihara:
		this->RR = new Tribology::AiharaR();
		break;
	case B4P_In::Tribology::RollingResistance::Fujiwara:
		this->RR = new Tribology::Fujiwara();
		break;
	case B4P_In::Tribology::RollingResistance::Houpert:
		this->RR = new Tribology::Houpert();
		break;
	case B4P_In::Tribology::RollingResistance::GoksemAihara:
		this->RR = new Tribology::GoksemAihara();
		break;
	}

	this->HZ = new Tribology::BrewHamrock();

	switch (FI.coulomb) {
	case B4P_In::Tribology::CoulombNothing:
		this->CL = new Tribology::CoulombNothing();
		break;
	case B4P_In::Tribology::Tangent:
		this->CL = new Tribology::Tangent();
		this->clmb_s = FI.coulomb_slope;
		break;
	}

	switch (FI.filmThickness) {
	case B4P_In::Tribology::FilmThicknessNothing:
		this->FT = new Tribology::FilmThicknessNothing();
		break;
	case B4P_In::Tribology::HamrockDowsonHc:
		this->FT = new Tribology::HamrockDowsonHc();
		break;
	case B4P_In::Tribology::HamrockDowsonHmin:
		this->FT = new Tribology::HamrockDowsonHmin();
		break;
	}

	this->TR = new Tribology::AiharaT();

	this->DF = new Tribology::Tsuji();

	switch(FI.hysteresis) {
	case B4P_In::Tribology::Kakuta:
		this->HY = new Tribology::Kakuta();
		this->fh = FI.hysteresis_factor;
		break;
	case B4P_In::Tribology::HysteresisNothing:
		this->HY = new Tribology::HysteresisNothing();
		this->fh = 0;
		break;
	}
	
	return;
}

// ���͂ɂ���ď����������Z�b�g����D
void B4P_BallRingPair::init_Lubrication(const B4P_In::Lubrication&FI) {

	this->lb_eta = FI.eta0;		// �S�x�iPa*s�j
	this->lb_beta = FI.beta0;	// ���x�S�x�W��
	this->lb_k = FI.k0;			// ���M�`����
	this->lb_alpha = FI.alpha0;	// ���͔S�x�W��
	this->lm = FI.lm0;			// ���j�X�J�X����

	return;
}

// ���͂ɂ���ĐڐG�ȉ~���̕ϐ����m�ۂ��Ă����D
void B4P_BallRingPair::init_Slice(int msmax) {

	for (int ig = 0; ig < 2; ig++) {
		int n = msmax;
		this->GV[ig].msmax = n;
		this->GV[ig].ps = new Vector3d[n];
		this->GV[ig].ratio_slice = new double[n];
		this->GV[ig].SL = new Slice[n];
		// �i�K�v�Ȃ��Ǝv�����C�j�X���C�X�Ђ̐��l��������
		this->init_Sliceparam(ig);
		// Hertz�ڐG���_����X���C�X�׏d�����߂�D
		Tribology::SliceForceRatio(n, this->GV[ig].ratio_slice);
	}
	return;
}

B4P_BallRingPair::B4P_BallRingPair() {
	for (int i = 0; i < 2; i++) {
		this->GV[i].ps = NULL;
		this->GV[i].ratio_slice = NULL;
	}
	return;
}

B4P_BallRingPair::~B4P_BallRingPair() {
	for (int i = 0; i < 2; i++) {
		if (this->GV[i].ps != NULL)
			delete[] this->GV[i].ps;
		if (this->GV[i].ratio_slice != NULL)
			delete[] this->GV[i].ratio_slice;
		if (this->GV[i].SL != NULL)
			delete[] this->GV[i].SL;
	}
	return;
}

// STEP1 �É�́@�����v�Z�F�ʍ��W�̑��
void B4P_BallOuterRingPair::set_eta_stf(const Vector3d & eta) {
	Vector3d x = this->RG->XZ_to_inecoord(eta);
	this->BL->x = x;
	return;
}


// STEP1 �É�́@�����v�Z�F�ʁ[���O�֊Ԃ̐ڐG�׏d�̌v�Z
void B4P_BallRingPair::get_Fstf(
	Vector3d&Fbi,	// out: [N]:	�]���̂����ցior�O�ցj����󂯂�� 
	Vector3d&Fib,	// out: [N]:	���ցior�O�ցj���]���̂���󂯂��
	Vector3d&Tib	// out:	[Nm]:	���ցior�O�ցj���]���̂���󂯂�g���N
) {
	Fbi = Vector3d::Zero();
	Fib = Tib = Vector3d::Zero();

	Vector3d bl_x = this->RG->to_mycoord(this->BL->x);
	Vector3d bl_v = this->RG->to_myvelocity(this->BL->v);
	Vector3d Tbi;

	// groove �� ball �̐ڐG�v�Z�D
	for (int i = 0; i < 2; i++) {

		// �K�v�Ȗ߂�l�̊m�ہD
		double dx, cos_alp, sin_alp, a, b, Pmax, lambda, fratio;
		dx = a = b = Pmax = lambda = fratio = 0;
		Vector3d er, eg, p, Fn, ur, us, Tr, Fs, Ts, Fi, Fb, Tb, Ti;
		er = eg = p = Fn = ur = us = Tr = Fs = Ts = Fi = Fb = Tb = Ti = Vector3d::Zero();

		// �ڂ��Ă��Ă��邩����D
		bool c = this->how_Contact(i, bl_x, er, eg, dx);
		if (c) {
			double Rx, Ry, k;
			this->calc_Hertz(i, bl_x, er, eg, dx, Rx, Ry, p, cos_alp, sin_alp, a, b, k);
			double F_norm = k * pow(dx, 1.5);
			Fb = F_norm * this->RG->to_inevector(-eg);
			Fi = F_norm * this->RG->to_inevector(eg);
			Pmax = 1.5 * F_norm / (Numeric::pi * a * b);
			Ti = this->RG->calc_Torque(p, Fi);	// �ʂɂ�����͂̔���p����a�ɂ�����g���N���v�Z����D
		}
		// �ڐG���Ă��Ȃ��ꍇ�ł��ア�΂˂Őڑ�
		else {
			double k = 1e-10;
			double F_norm = k * dx;
			Fb = F_norm * this->RG->to_inevector(-eg);
			Fi = F_norm * this->RG->to_inevector(eg);
			Ti = this->RG->calc_Torque(p, Fi);	// �ʂɂ�����͂̔���p����a�ɂ�����g���N���v�Z����D
		}
		// �v�Z���ʂ������o�ϐ��ɕۑ�����D
		//this->save(c, i, cos_alp, sin_alp, a, b, p, Fn, Pmax, us, ur, Tr, lambda, fratio, Fs, Ts, Fb, Tb, Ti);

		// �a0�ƍa1�̍��v���o�͂���D
		Fbi += Fb;		// Fbi : �{�[��(b)�������O(i)����󂯂�́D
		Fib += Fi;		// Fib : �����O(i)���{�[��(b)����󂯂�́D
		Tib += Ti;
	}

	return;
}


// �e�X���C�X�̉׏d���o��
//VectorXd B4P_BallRingPair::make_SliceParam(void) {
//	int NParam = 15;
//	int Nslice = this->GV[0].msmax + this->GV[1].msmax;
//	VectorXd Param(Nslice * NParam);
//	for (int i = 0; i < 2; i++) {
//		int Ni = i * NParam * this->GV[0].msmax;
//		for (int k = 0; k < this->GV[i].msmax; k++) {
//			int Nik = Ni + NParam * k;
//			for (int j = 0; j < 3; j++) {
//				int Nikj = Nik + j;
//				Param(Nikj + 0) = this->GV[i].ps[k][j];
//				Param(Nikj + 3) = this->GV[i].Fs_slice[k][j];
//				Param(Nikj + 6) = this->GV[i].Ts_slice[k][j];
//				Param(Nikj + 9) = this->GV[i].us_slice[k][j];
//			}
//			Param(Nik + 12) = this->GV[i].fn_slice[k];
//			Param(Nik + 13) = this->GV[i].mutr_slice[k];
//			Param(Nik + 14) = this->GV[i].mucl_slice[k];
//		}
//	}
//	return Param;
//}



//Vector3d ub = this->BL->surface_velocity(p);
//Vector3d ui = this->RG->surface_velocity(p);
//Vector3d ur = 0.5 * (ub + ui);


// ���ނ̍ޗ��ƐڐG����������C�_���s���O�W�������߂�D"15_����`�o�l�ƏՓ˂̃��f����"�Q�l�D
//// �v�Z����"���{�f���Ȋw�, 2004, 13.4: 233-240."�̔���`�΂˂ւ̎�v�Z�ɂ��K���D 
//double B4P_BallRingPair::get_c
//(				// out:	[Ns/m]:	�_���s���O�W���D
//	double k	// in:	[N/m]:	�ڐG�������D
//) {
//	double c = 2 * this->zeta * sqrt(1.5 * this->m * k);
//	return c;
//}




//// ���S�͂��Z�o���郁�\�b�h�i�É�͗p�j�߂�l�F�������W�n[N]
//Vector3d B4P_BallRingPair::get_centf() {
//
//	// ���O�֍��W�n(�����O�֒��S�)�̋ʂ̈ʒu�x�N�g��
//	Vector3d r_rg = this->RG->to_mycoord(this->BL->x);
//
//	// x�����͊֌W�Ȃ��̂�0�ɐݒ�
//	r_rg.x() = 0;
//
//	// r_�̃m�����ƒP�ʃx�N�g�����`
//	double rrg_nor = r_rg.norm();
//	Vector3d e_rg = r_rg / rrg_nor;
//
//	// ���O�֍��W�n�̋ʑ��x�̐ڐ������������Z�o
//	Vector3d v_rg = this->RG->to_mycoord(this->BL->v);
//	v_rg.x() = 0;
//	double vrg_th = (v_rg.cross(e_rg)).norm();
//
//	// ���S�͂̌v�Z��(F=m*v^2/r)��K�p���C���O�֍��W�n��̉��S�͂��o��
//	Vector3d F_rg = e_rg * this->BL->m * vrg_th * vrg_th / rrg_nor;
//
//	// ���O�֍��W�n�̉��S�͂��������W�n�ɕϊ�
//	Vector3d F_inr = this->RG->to_inevector(F_rg);
//
//	// �������W�n�ɂ����鉓�S�͂��o��
//	return F_inr;
//}


//// �]���葬�x���Z�o���郁�\�b�h�i�������W�n�j(��xx�āC�����Ԉ���Ă���)
//void B4P_BallRingPair::get_us_ur2(
//	const Vector3d&p,		// in:	[-]:	�ڐG�_�ʒu�D�i�������W�n�j
//	const Vector3d&dpdt,		// in:	[-]:	�ڐG�_�ʒu�D�i�������W�n�j
//	Vector3d & us,			// out:	[-]:	���葬�x�D�i�������W�n�j
//	Vector3d & ur			// out:	[-]:	�]���葬�x�D�i�������W�n�j
//) {
//	Vector3d dx = p - this->BL->x;
//	Vector3d e = dx.normalized();
//	Vector3d rwb = this->BL->w.cross(dx);
//	Vector3d rwr = this->RG->w.cross(p - this->RG->x);
//	Vector3d bl_v = this->BL->v - this->RG->v;
//	//Vector3d e_a = this->RG->get_ax();
//
//	Vector3d RB = this->BL->x - this->RG->x;
//	
//
//
//	Vector3d ur_ = 0.5 * (this->BL->v + this->RG->v - 2 * dpdt + rwb + rwr);
//	ur = ur_-e * e.dot(ur_);
//
//	Vector3d ub = this->BL->v + rwb;
//	Vector3d ui = this->RG->v + rwr;
//	Vector3d us_ = ub - ui;
//	us = us_;// -e //* e.dot(us_);
//
//	return;
//}