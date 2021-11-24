/*******************************************************************************
!								"BS_BallSpiralPair.cpp"
!													2020/01/08	[Core-T]	���
!	�{�[���Ɖ~�����ށi�i�b�g�C�V���t�g�j�̃y�A�N���X�D
!	���C�W���̂悤�ȕ��ފԕs�ϗʂ̕ێ��̑��C�ڐG�׏d�Ȃǂ̌v�Z���s���D
!
!*******************************************************************************/
#include "BS_BallCylinderPair.h"

// �ǂ̉~���̂ǂ̗����ɕR�Â������`����D
void BS_BallCylinderPair::link(Ball*BL, BS_Cylinder*CY, int is) {
	this->BL = BL;
	this->CY = CY;
	this->iSP = is;

	return;
}

// ���݂̏�Ԃł̋ʂ̐ڐG���肨��сC���̂Ƃ��̐ڐG��Ԃ��v�Z����D
bool BS_BallCylinderPair::how_Contact(
	int i,					// in:	[-]		�ڐG�a�ԍ��D0 or 1�D
	const Vector2d&BL_eta,	// in:	[m]		�a���s�f�ʂɂ�����{�[���̍��W�D�i���_�FPCD�������S�C��-�č��W�j
	Vector2d&e,				// out:	[-]:	�a���S���猩���{�[�������x�N�g���D�{�[���̎󂯂�ڐG�׏d�����̋t�����D(��-�č��W)
	double&dx				// out:	[-]		�e���ڋߗʁD���ŐڐG�D
) {
	// �a���S�_����{�[�����S�Ɍ������x�N�g���̎Z�o�D
	Vector2d gb = BL_eta - this->CY->SP[iSP].GV[i].eta;
	double gb_norm = gb.norm();

	// �􉽊w����{�[���H�����ݗʂ̎Z�o�D
	dx = gb_norm - this->CY->SP[iSP].GV[i].r + this->BL->r;

	// �a���S�_�ƃ{�[�����S����v���Ă�����ڂ��Ă��Ȃ����Ƃɂ���D�i0���Z�h�~�j
	if (gb_norm == 0.0)
		return false;

	e = gb / gb_norm;

	// �{�[�����S���g�[���X����o�Ă�����ڂ��Ă��Ȃ����Ƃɂ���D
	if (gb_norm > this->CY->SP[iSP].GV[i].r) {
		return false;
	}

	// �H�����ݗʂ����l�̏ꍇ�ڐG���Ă��Ȃ��D
	if (dx <= 0)
		return false;

	return true;
}

// �ڐG�_�������߂郁�\�b�h�D
int BS_BallCylinderPair::num_Contact() {

	Vector3d th_eta = this->CY->to_etacoord(iSP, this->BL->x);
	Vector2d eta = Vector2d(th_eta[1], th_eta[2]);		// �ʂ̃�-�č��W�i�������W�j

	Vector2d e;
	double dx;

	int nc = 0;
	for (int i = 0; i < 2; i++)
		nc += int(this->how_Contact(i, eta, e, dx));

	return nc;
}

// Hertz�̌v�Z����C�׏d��Ԃ����\�b�h�D�_���p�ɂ�錸���Ȃ��D
double BS_BallCylinderPair::calc_Hertz(
	int ig,					// in:	[-]		�ڐG�a�ԍ��D0 or 1�D
	const Vector2d&bl_eta,	// in:	[m]		�������W�n�ł̃{�[���ʒu�D
	const Vector2d&e,		// in:	[-]:	�a���S���猩���{�[�������x�N�g���D�{�[���̎󂯂�ڐG�׏d�����̋t�����D
	double dx,				// in:	[-]		�e���ڋߗʁD���ŐڐG�D
	double&Rx,				// out:	[m]:	�]��������ȗ����a�D
	double&Ry,				// out:	[m]:	�������ȗ����a�D
	Vector2d&p,				// out:	[m]:	�ڐG�_�ʒu�i�a���s�n�j�D
	double&cos_alp,			// out:	[-]:	�ڐG�p��cos�D
	double&sin_alp,			// out:	[-]:	�ڐG�p��sin�D
	double&a,				// out:	[m]:	�ڐG�ȉ~���a�D
	double&b,				// out:	[m]:	�ڐG�ȉ~�Z�a�D
	double&k,				// out:	[N/x^1.5]:	����`�����D
	Vector2d&rho			// out:	[1/m]:	�˂����̋ȗ�
) {
	cos_alp = e[0];
	sin_alp = e[1];

	rho = this->CY->SP[this->iSP].get_rho(cos_alp, ig);

	Rx = 1.0 / (this->BL->r_inv + rho[0]);
	Ry = 1.0 / (this->BL->r_inv + rho[1]);

	// BrewHamrock �̋ߎ������獄���C�ڐG�ȉ~�����߂�D 
	this->TB.HZ->calc(Rx, Ry, dx, this->E, k, a, b);
	double F_norm_ = k * pow(dx, 1.5);

	// �����ȗ����a����ڐG�ȉ~���S�����߂�D
	double r0 = Numeric::EffectiveCenter(this->BL->r, -1 / rho[1], dx);
	p = bl_eta + r0 * e;

	return F_norm_;
}

// Hertz�̌v�Z����C�׏d��Ԃ����\�b�h�D�_���p�ɂ�錸������D
double BS_BallCylinderPair::calc_DynamicHertz(
	int ig,					// in:	[-]		�ڐG�a�ԍ��D0 or 1�D
	const Vector2d&bl_eta,	// in:	[m]		�������W�n�ł̃{�[���ʒu�D
	const Vector2d&bl_etav,	// in:	[m]		�������W�n�ł̃{�[�����x�D
	const Vector2d&e,		// in:	[-]:	�a���S���猩���{�[�������x�N�g���D�{�[���̎󂯂�ڐG�׏d�����̋t�����D
	double dx,				// in:	[-]		�e���ڋߗʁD���ŐڐG�D
	double&Rx,				// out:	[m]:	�]��������ȗ����a�D
	double&Ry,				// out:	[m]:	�������ȗ����a�D
	Vector2d&p,				// out:	[m]:	�ڐG�_�ʒu�i�������W�n�j�D
	double&cos_alp,			// out:	[-]:	�ڐG�p��cos�D
	double&sin_alp,			// out:	[-]:	�ڐG�p��sin�D
	double&a,				// out:	[m]:	�ڐG�ȉ~���a�D
	double&b,				// out:	[m]:	�ڐG�ȉ~�Z�a�D
	Vector2d&rho			// out:	[1/m]:	�˂����̋ȗ�
) {
	// BrewHamrock �̋ߎ������獄���C�ڐG�ȉ~�����߂�D 
	double k;
	double F_norm_ = this->calc_Hertz(ig, bl_eta, e, dx, Rx, Ry, p, cos_alp, sin_alp, a, b, k, rho);

	// ����`��������C����`�_���p�W�� c ���Z�o�D
	double c = 2 * this->GV[ig].zeta * sqrt(1.5 * this->m * k);

	// �{�[�����ǖʂɌ������Ă������x�������Z�o�D
	double v = bl_etav.dot(e);

	// �����͂��׏d�ɉ��Z�D
	double Ftotal = F_norm_ + c * v * pow(dx, 0.25);

	// �׏d�����l�����Ȃ��悤�C�������l�ŉ�����ݒ�D
	return std::max(Ftotal, 1e-20);
}

// ���݂̏�Ԃ���C�{�[���Ɖ~���̐ڐG�׏d�݂̂����߂郁�\�b�h�D
void BS_BallCylinderPair::get_F0(
	Vector2d&Fbc,	// out:	[N]:	�{�[��(b)���~��(c)����󂯂�́i�����������W�C0:�Ő����C1:�Đ����j�D
	Vector3d&Fcb,	// out:	[N]:	�~��(c)���{�[��(b)����󂯂�́i�������W�n�j�D
	Vector3d&Tcb	// out:	[m]:	�~��(c)���{�[��(b)����󂯂�g���N�i�������W�n�j�D
) {
	Vector3d th_eta = this->CY->to_etacoord(iSP, this->BL->x);
	double th = th_eta[0];								// �ʂ̈ʑ��p[rad]
	Vector2d eta = Vector2d(th_eta[1], th_eta[2]);		// �ʂ̃�-�č��W�i�������W�j
	Matrix3d xyz2eta = this->CY->get_xyz2eta(iSP, th);

	Fbc = Vector2d::Zero();
	Fcb = Tcb = Vector3d::Zero();

	for (int ig = 0; ig < 2; ig++) {

		Vector3d p_, Fcb_;
		p_ = Fcb_ = Vector3d::Zero();
		double F_norm, a;
		F_norm = a = 0.0;

		double dx;
		Vector2d e, Fbc_;
		Fbc_ = Vector2d::Zero();

		bool is_contact = this->how_Contact(ig, eta, e, dx);

		if (is_contact) {
			double Rx, Ry, cos_alp, sin_alp, b, k;
			Vector2d p, rho;
			F_norm = this->calc_Hertz(ig, eta, e, dx, Rx, Ry, p, cos_alp, sin_alp, a, b, k, rho);
			Fbc_ = -e * F_norm;
			Fcb_ = this->CY->to_inertialvector(Vector3d(0.0, -Fbc_[0], -Fbc_[1]), xyz2eta);
			p_ = this->CY->to_inertialcoord(iSP, Vector3d(th, p[0], p[1]));
		}
		this->save_F0(ig, p_, F_norm, th_eta, a, Fbc_);
		Fbc += Fbc_;
		Fcb += Fcb_;
		Tcb += this->CY->calc_Torque(p_, Fcb_);
	}
	return;
}

void BS_BallCylinderPair::save_F0(int ig, const Vector3d&p, double F, const Vector3d&eta, double a, const Vector2d&Fbc) {

	this->SV.eta = eta;

	this->SV.GV[ig].p = p;
	this->SV.GV[ig].F = F;
	this->SV.GV[ig].a = a;
	this->SV.GV[ig].Fbc = Fbc;

	return;
}

// ���݂̏�Ԃ���C�{�[���Ɖ~���̐ڐG�׏d�݂̂����߂郁�\�b�h�D
void BS_BallCylinderPair::get_F1(
	bool direction,		// in:	[-]		�~�����猩���{�[���i�s�����x�N�g��(+x:false�C-x:true)
	Vector2d&Fbc,		// out:	[N]:	�{�[��(b)���~��(c)����󂯂�́i�����������W�C�������E�a�����j�D
	Vector3d&Fcb,		// out:	[N]:	�~��(c)���{�[��(b)����󂯂�́i�������W�n�j�D
	Vector3d&Tcb		// out:	[Nm]:	�~��(c)���{�[��(b)����󂯂�g���N�i�������W�n�j�D
) {
	Vector3d e = this->get_e();

	Vector3d th_eta = this->CY->to_etacoord(iSP, this->BL->x);
	double th = th_eta[0];
	Vector2d eta = Vector2d(th_eta[1], th_eta[2]);
	Matrix3d xyz2eta = this->CY->get_xyz2eta(iSP, th);

	Fbc = Vector2d::Zero();
	Fcb = Tcb = Vector3d::Zero();

	for (int ig = 0; ig < 2; ig++) {

		Vector3d p_, Fcb_;
		p_ = Fcb_ = Vector3d::Zero();

		double F_norm, a, dx;
		F_norm = a = 0.0;

		Vector2d e_, Fbc_;
		Fbc_ = Vector2d::Zero();

		bool is_contact = this->how_Contact(ig, eta, e_, dx);

		if (is_contact) {
			double Rx, Ry, cos_alp, sin_alp, b, k;
			Vector2d p, rho;
			F_norm = this->calc_Hertz(ig, eta, e_, dx, Rx, Ry, p, cos_alp, sin_alp, a, b, k, rho);
			Fbc_ = -e_ * F_norm;

			// �i�b�g/�V���t�g���猩���ʐi�s�����̌������疀�C�͂̕��������ߑł��Ōv�Z
			Vector2d muFb;
			double mu = this->GV[ig].mu;

			if (direction)
				muFb = Vector2d(-mu * Fbc_[1], mu * Fbc_[0]);
			else
				muFb = Vector2d(mu * Fbc_[1], -mu * Fbc_[0]);

			Fbc_ += muFb;
			p_ = this->CY->to_inertialcoord(iSP, Vector3d(th, p[0], p[1]));

			//// �Ă��Ɓ�
			//double k_ = -1e6;
			//Fbc_ = -e_ * k_ * dx;
			//// �Ă��Ɓ�
		}
		//else {
		//	// �Ă��Ɓ�
		//	double k = -1e6;
		//	Fbc_ = -e_ * k * dx;
		//	// �Ă��Ɓ�
		//}
		Fcb_ = this->CY->to_inertialvector(Vector3d(0.0, -Fbc_[0], -Fbc_[1]), xyz2eta);

		this->save_F0(ig, p_, F_norm, th_eta, a, Fbc_);
		Fbc += Fbc_;
		Fcb += Fcb_;
		Tcb += this->CY->calc_Torque(p_, Fcb_);
	}
	return;
}

// ���x�ˑ��̂Ȃ����C�����߂郁�\�b�h�D
void BS_BallCylinderPair::get_F2(
	Vector3d&vF,		// out:	[Nm/s]	�i�b�g�i�s�������薀�C�́D
	Vector3d&vT			// out:	[Nm2/s]	���薀�C�ɂ��g���N�i�������W�n�j�D
) {
	vF = vT = Vector3d::Zero();

	for (int ig = 0; ig < 2; ig++) {

		Vector3d p = this->SV.GV[ig].p;
		Vector3d v = this->get_us(p);
		double   F = this->SV.GV[ig].F;
		Vector3d vF_ = v * F;
		Vector3d vT_ = this->BL->calc_Torque(p, vF_);
		vF += vF_;
		vT += vT_;
	}
	return;
}

// ���݂̏�Ԃ���C�{�[���Ɖ~���̐ڐG�׏d�݂̂����߂郁�\�b�h�D
void BS_BallCylinderPair::get_dyn_F0(
	bool direction,		// in:	[-]		�~�����猩���{�[���i�s�����x�N�g��(+x:false�C-x:true)
	Vector3d&Fbc,		// out:	[N]:	�{�[��(b)���~��(c)����󂯂�́i�������W�n�j�D
	Vector3d&Fcb,		// out:	[N]:	�~��(c)���{�[��(b)����󂯂�́i�������W�n�j�D
	Vector3d&Tcb		// out:	[Nm]:	�~��(c)���{�[��(b)����󂯂�g���N�i�������W�n�j�D
) {
	Vector3d th_eta = this->CY->to_etacoord(iSP, this->BL->x);
	double th = th_eta[0];
	Vector2d eta = Vector2d(th_eta[1], th_eta[2]);
	Matrix3d xyz2eta = this->CY->get_xyz2eta(iSP, th);

	Fbc = Fcb = Tcb = Vector3d::Zero();

	for (int ig = 0; ig < 2; ig++) {

		Vector3d p_, Fcb_;
		p_ = Fcb_ = Vector3d::Zero();

		double F_norm, a, dx;
		F_norm = a = 0.0;

		Vector2d e_, Fbc_;
		Fbc_ = Vector2d::Zero();

		bool is_contact = this->how_Contact(ig, eta, e_, dx);

		if (is_contact) {
			double Rx, Ry, cos_alp, sin_alp, b, k;
			Vector2d p, rho;
			F_norm = this->calc_Hertz(ig, eta, e_, dx, Rx, Ry, p, cos_alp, sin_alp, a, b, k, rho);
			Fbc_ = -e_ * F_norm;

			// �i�b�g/�V���t�g���猩���ʐi�s�����̌������疀�C�͂̕��������ߑł��Ōv�Z
			Vector2d muFb;
			double mu = this->GV[ig].mu;

			if (direction)
				muFb = Vector2d(-mu * Fbc_[1], mu * Fbc_[0]);
			else
				muFb = Vector2d(mu * Fbc_[1], -mu * Fbc_[0]);

			Fbc_ += muFb;
			p_ = this->CY->to_inertialcoord(iSP, Vector3d(th, p[0], p[1]));
		}
		Fcb_ = this->CY->to_inertialvector(Vector3d(0.0, -Fbc_[0], -Fbc_[1]), xyz2eta);

		Fbc += Fbc_;
		Fcb += Fcb_;
		Tcb += this->CY->calc_Torque(p_, Fcb_);
	}
	return;
}

// Hertz�ڐG�����C���܂߂��S�Ă̗͂�Ԃ��D��ɓ���͗p�D
void BS_BallCylinderPair::get_FT(
	Vector3d&Fbc,	// out
	Vector3d&Tbc,	// out
	Vector3d&Tcb,	// out
	Vector3d&Fs,	// out	�~�����ނɂ����銊��ɂ��׏d�D
	Vector3d&Ts		// out	�~�����ނɂ����銊��ɂ��g���N�D
) {
	Fbc = Tbc = Tcb = Fs = Ts = Vector3d::Zero();

	Vector3d th_eta = this->CY->to_etacoord(this->iSP, this->BL->x);
	Vector2d eta = Vector2d(th_eta[1], th_eta[2]);
	double th = th_eta[0];
	Matrix3d xyz2eta = this->CY->get_xyz2eta(this->iSP, th);

	Vector3d BL_v = this->BL->v;
	Vector3d etav = this->CY->to_etavelocity(BL_v, xyz2eta);
	Vector2d etav_ = Vector2d(etav[1], etav[2]);

	Vector3d exai = this->CY->to_inertialvector(Vector3d::UnitX(), xyz2eta);	// �������W�n�̍a�f�ʐi�s����

	for (int ig = 0; ig < 2; ig++) {
		double dx, Rx, Ry, cos_alp, sin_alp, a, b, F_norm, Pmax, lambda, fratio, h;
		dx = Rx = Ry = cos_alp = sin_alp = a = b = F_norm = Pmax = lambda = fratio = h = 0.0;
		Vector2d e, p, rho;
		e = p = rho = Vector2d::Zero();
		Vector3d p_, Fn, us, ur, Tbr, Fs_, Ts_, Fb, Tb, Ti, Tcs_;
		p_ = Fn = us = ur = Tbr = Fs_ = Ts_ = Fb = Tb = Ti = Tcs_ = Vector3d::Zero();

		bool is_contact = this->how_Contact(ig, eta, e, dx);

		if (is_contact) {
			F_norm = this->calc_DynamicHertz(iSP, eta, etav_, e, dx, Rx, Ry, p, cos_alp, sin_alp, a, b, rho);
			p_ = this->CY->to_inertialcoord(iSP, Vector3d(th, p[0], p[1]));
			Fn = this->CY->to_inertialvector(
				Vector3d(0.0, F_norm * -e[0], F_norm * -e[1]),
				xyz2eta);

			Pmax = 1.5 * F_norm / (Numeric::pi * a * b);
			ur = this->get_ur(p_);
			us = this->get_us(p_);
			double ur_norm = ur.norm();

			// �{�[���ɂ����銊�薀�C�R�͂̎Z�o�D
			h = this->TB.FT->calc(F_norm, this->E, Rx, Ry, this->LB.alpha, this->LB.eta, ur_norm, this->LB.lm, b);
			h *= Tribology::ErtelGrubin(this->LB.eta, this->LB.beta, this->LB.k, ur_norm);
			lambda = h / this->GV[ig].sigma;
			fratio = Tribology::ForceRatio(lambda);
			this->calc_Sliding(ig, th, p_, exai, a, b, Pmax, ur_norm, F_norm, fratio, rho, Fs_, Ts_, Tcs_);

			// �{�[���ɂ�����]���薀�C�R�͂̎Z�o�D�i�������W�n�j
			double Trr_norm = this->TB.RR->calc(Rx, Ry, this->BL->D, a, b, F_norm, this->E, ur_norm, fratio, this->LB.eta, this->LB.alpha, this->LB.beta, this->LB.k, this->LB.lm);
			double Trh_norm = this->TB.HY->calc(this->TB.fh, b, F_norm);
			Tbr = (Trr_norm + Trh_norm) * this->BL->calc_TorqueDirection(p_, -ur);

			Fb = Fn + Fs_;
			Tb = Tbr + Ts_;
			Ti = this->CY->calc_Torque(p_, -Fn) + Tcs_ - Tbr;
			Fs -= Fs_;
			Ts += Tcs_ - Tbr;
		}
		else {
			// �ڐG���Ȃ��Ƃ��C�X���C�X�Ђ̌��ʂ����ׂ�0�ɂ���
			this->init_Sliceparam(ig);
		}
		this->save(ig, th_eta, Fn, Fs_, us, ur, dx, cos_alp, sin_alp, a, b, h, p_, lambda, fratio, Pmax);

		// �a0�ƍa1�̍��v���o�͂���D
		Fbc += Fb;		// Fbc : �{�[��(b)���~��(c)����󂯂�́D
		Tbc += Tb;		// Tbc : �{�[��(b)���~��(c)����󂯂郂�[�����g�D
		Tcb += Ti;		// Tcb : �~��(c)���{�[��(b)����󂯂郂�[�����g�D
	}
	return;
}

// �]���葬�x���Z�o���郁�\�b�h�i�������W�n�j
Vector3d BS_BallCylinderPair::get_ur(
	const Vector3d&p		// �ڐG�_�ʒu�D�i�������W�n�j
) {
	Vector3d rwb = this->BL->w.cross(p - this->BL->x);
	Vector3d rwc = this->CY->w.cross(p - this->CY->x);
	Vector3d ur = 0.5 * (3 * this->CY->v - this->BL->v + rwb + rwc);
	return ur;
}

// ���葬�x���Z�o���郁�\�b�h�i�������W�n�j
Vector3d BS_BallCylinderPair::get_us(
	const Vector3d&p		// �ڐG�_�ʒu�D�i�������W�n�j
) {
	Vector3d ub = this->BL->surface_velocity(p);
	Vector3d ui = this->CY->surface_velocity(p);
	Vector3d us = ub - ui;
	return us;
}

// 
Vector3d BS_BallCylinderPair::get_eta(void) {
	Vector3d eta = this->CY->to_etacoord(iSP, this->BL->x);
	return eta;
}

// ���薀�C�����߂郁�\�b�h�D
void BS_BallCylinderPair::calc_Sliding
(
	int ig,				// in:	[-]		�ڐG�a�ԍ��D0 or 1�D
	double th,			// in:	[rad]:	�����ʑ��p�D
	const Vector3d&p,	// in:	[m]:	�ڐG�_�ʒu�i�������W�n�j�D
	const Vector3d&xai,	// in:	[m]:	�a�f�ʂɂ����闆���i�s�����x�N�g���D�i�������W�n�j
	double a,			// in:	[m]:	�ڐG�ȉ~���a�D
	double b,			// in:	[m]:	�ڐG�ȉ~�Z�a�D
	double Pmax,		// in:	
	double ur_norm,		// in:	
	double F_norm,		// in:	
	double fratio,		// in:	
	const Vector2d&rho,	// in:	[1/m]:	�˂����̋ȗ��D
	Vector3d&Fs,		// out:	[N]:	�{�[���ɂ����銊�薀�C�́i�������W�n�j
	Vector3d&Ts,		// out:	[Nm]:	�{�[���ɂ����銊�薀�C�g���N�i�������W�n�j
	Vector3d&Tcs		// out:	[Nm]:	�a�ɂ����銊�薀�C�g���N�i�������W�n�j
) {
	// �ȗ��������a����ڐG�ȉ~�ʒu�̃��b�V�������߂�D
	Fs = Ts = Tcs = Vector3d::Zero();
	double rhox = this->BL->r_inv - rho[1];
	double rhoy = this->BL->r_inv - rho[0];
	Vector2d rhom(rhox, rhoy);

	this->CY->calc_slice(iSP, ig, th, p, xai, a, b, this->GV[ig].xy, this->GV[ig].n, rhom, this->GV[ig].ps);

	// �ڐG�ȉ~���X���C�X���Ă��ꂼ��̊��薀�C�͂�ώZ����D
	for (int j = 0; j < this->GV[ig].n; j++) {
		Vector3d ps = this->GV[ig].ps[j];
		Vector3d us_ = this->get_us(ps);
		double us_norm = us_.norm();
		double mu_tr = this->TB.TR->calc(this->LB.eta, Pmax, us_norm, ur_norm);
		double mu_cl = this->TB.CL->calc(this->GV[ig].mu, us_norm, ur_norm, this->TB.cs);
		double f_arr = this->GV[ig].r[j] * F_norm;
		double Fs_norm = mu_tr * fratio * f_arr + mu_cl * (1.0 - fratio) * f_arr;
		Vector3d Fs_ = Fs_norm * -us_ / us_norm;
		Vector3d Ts_ = this->BL->calc_Torque(ps, Fs_);
		Vector3d Tcs_ = this->CY->calc_Torque(ps, -Fs_);
		Fs += Fs_;
		Ts += Ts_;
		Tcs += Tcs_;
		this->save_Slice(ig, j, f_arr, Fs_, Ts_, mu_cl, mu_tr, us_, this->GV[ig].ps[j]);
	}
	return;
}

// �X���C�X�v�Z���ʂ������o�ϐ��Ɋi�[
void BS_BallCylinderPair::save_Slice(
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
	this->SV.GV[i].SL[j].f_arr = f_arr;
	this->SV.GV[i].SL[j].fs = Fs_;			// ���薀�C��[N](�������W�n)
	this->SV.GV[i].SL[j].ts = Ts_;			// ���薀�C�g���N[Nm](�������W�n)
	this->SV.GV[i].SL[j].mu_cl = mu_cl;		// �N�[�������C�W��[-]
	this->SV.GV[i].SL[j].mu_tr = mu_tr;		// �g���N�V�����W��[-]
	this->SV.GV[i].SL[j].us = us_;			// ���葬�x[m/s](�������W�n)
	this->SV.GV[i].SL[j].ps = ps_;			// �X���C�X�Ђ̑��Έʒu[m](�������W�n)
	return;
}

void BS_BallCylinderPair::init(
	const BS_In::BallCylinderPair & BCP,
	const BS_In::Tribology & tribology,
	const BS_In::Oil & oil
) {
	for (int k = 0; k < 2; k++) {

		switch (tribology.ellipse) {
		case BS_In::Tribology::Ellipse::Mesh1d: {
			int n0 = tribology.ellipse_mesh[0];
			this->GV[k].n = tribology.ellipse_mesh[0];
			this->GV[k].r = new double[this->GV[k].n];
			this->GV[k].xy = new Vector2d[this->GV[k].n];
			Tribology::SliceForceRatio(this->GV[k].n, this->GV[k].r);
			for (int i = 0; i < n0; i++)
				this->GV[k].xy[i] = Vector2d(double(2 * i - n0 + 1) / n0, 0);

			break;
		}
		case BS_In::Tribology::Ellipse::Mesh2d: {
			int n0 = tribology.ellipse_mesh[0];
			int n1 = tribology.ellipse_mesh[1];
			double*r_temp = new double[n0];
			Tribology::SliceForceRatio2d(n0, n1, r_temp);
			this->GV[k].n = n0 * n1;
			this->GV[k].r = new double[this->GV[k].n];
			this->GV[k].xy = new Vector2d[this->GV[k].n];
			for (int i = 0; i < n0; i++) {
				double r0 = (i + 0.5) / n0;
				for (int j = 0; j < n1; j++) {
					int ij = i * n1 + j;
					double th = (j + 0.5) / n1 * 2 * Numeric::pi;

					this->GV[k].r[ij] = r_temp[i];
					this->GV[k].xy[ij] = r0 * Vector2d(cos(th), sin(th));
				}
			}
			break;
		}
		}
		this->GV[k].ps = new Vector3d[this->GV[k].n];

		this->GV[k].sigma = Vector2d(
			this->BL->sigmap,
			this->CY->SP[iSP].GV[k].sigma
		).norm();		//�����e��

		this->GV[k].mu = BCP.groove[k].mu;
		this->GV[k].zeta = BCP.groove[k].zeta;

		this->SV.GV[k].SL = new Save::Groove::Slice[this->GV[k].n];

		this->init_Sliceparam(k);
	}
	this->m = this->BL->m; // 1.0 / (1.0 / this->BL->m  + 1.0 / this->CY->m);//���Z����
	this->E = Tribology::ReducedYoung(this->BL->E, this->BL->nu, this->CY->E, this->CY->nu);	//���������O��

	// �g�p����g���C�{���W�[���̏������D
	this->init_Tribology(tribology);
	this->init_Oil(oil);

	return;
}

// �X���C�X�Ђ̌��ʂ�0���i�[
void BS_BallCylinderPair::init_Sliceparam(int i) {
	for (int j = 0; j < this->GV[i].n; j++) {
		this->SV.GV[i].SL[j].f_arr = 0;
		this->SV.GV[i].SL[j].fs = Vector3d::Zero();	// ���薀�C��[N](�������W�n)
		this->SV.GV[i].SL[j].ts = Vector3d::Zero();	// ���薀�C�g���N[Nm](�������W�n)
		this->SV.GV[i].SL[j].mu_cl = 0;				// �N�[�������C�W��[-]
		this->SV.GV[i].SL[j].mu_tr = 0;				// �g���N�V�����W��[-]
		this->SV.GV[i].SL[j].us = Vector3d::Zero();	// ���葬�x[m/s](�������W�n)
		this->SV.GV[i].SL[j].ps = Vector3d::Zero();	// ���葬�x[m/s](�������W�n)			
	}
}

Vector3d BS_BallCylinderPair::get_etav0(void) {

	Vector3d th_eta = this->CY->to_etacoord(this->iSP, this->BL->x);
	Matrix3d xyz2eta = this->CY->get_xyz2eta(this->iSP, th_eta[0]);

	Vector3d etav = this->CY->to_etavelocity(this->BL->v, xyz2eta);

	return etav;
}

Vector3d BS_BallNutPair::get_eta0(void) {

	return this->CY->to_etacoord(this->iSP, this->BL->x);
}

void BS_BallNutPair::set_eta0(const Vector2d & eta) {

	Vector3d eta_ = Vector3d(this->th0, eta[0], eta[1]);
	Vector3d x = this->CY->to_inertialcoord(iSP, eta_);
	this->BL->x = x;

	return;
}

// �i�b�g�̍a���s�f�ʂ́C2�������������o�����\�b�h�D
Vector2d BS_BallNutPair::get_etavector0(const Vector3d & x) {

	Matrix3d xyz2eta = this->CY->get_xyz2eta(iSP, th0);
	Vector3d eta = this->CY->to_etavector(x, xyz2eta);
	Vector2d eta_ = Vector2d(eta[1], eta[2]);

	return eta_;
}


void BS_BallShaftPair::set_etap(const Vector2d & eta) {

	Vector3d eta_ = Vector3d(this->thp, eta[0], eta[1]);
	Vector3d x = this->CY->to_inertialcoord(iSP, eta_);
	this->BL->x = x;

	return;
}

// �i�b�g�̍a���s�f�ʂ́C2�������������o�����\�b�h�D
Vector2d BS_BallShaftPair::get_etavectorp(const Vector3d & x) {

	Matrix3d xyz2eta = this->CY->get_xyz2eta(iSP, thp);
	Vector3d eta = this->CY->to_inertialvector(x, xyz2eta);
	Vector2d eta_ = Vector2d(eta[1], eta[2]);

	return eta_;
}

void BS_BallCylinderPair::init_Tribology(const BS_In::Tribology & tribology) {

	// �g�p����g���C�{���W�[���̏������D
	this->TB.HZ = new Tribology::BrewHamrock();

	switch (tribology.rollingresistance) {
	case BS_In::Tribology::RollingResistance::RollingResistanceNothing:
		this->TB.RR = new Tribology::RollingResistanceNothing();
		break;
	case BS_In::Tribology::RollingResistance::Aihara:
		this->TB.RR = new Tribology::AiharaR();
		break;
	case BS_In::Tribology::RollingResistance::Fijiwara:
		this->TB.RR = new Tribology::Fujiwara();
		break;
	case BS_In::Tribology::RollingResistance::Houpert:
		this->TB.RR = new Tribology::Houpert();
		break;
	case BS_In::Tribology::RollingResistance::GoksemAihara:
		this->TB.RR = new Tribology::GoksemAihara();
		break;
	}

	switch (tribology.coulomb) {
	case BS_In::Tribology::Coulomb::CoulombNothing:
		this->TB.CL = new Tribology::CoulombNothing();
		break;
	case BS_In::Tribology::Coulomb::Tangent:
		this->TB.CL = new Tribology::Tangent();
		this->TB.cs = tribology.coulomb_slope;
		break;
	case BS_In::Tribology::Coulomb::SimpleCoulomb:
		this->TB.CL = new Tribology::SimpleCoulomb();
		break;
	}

	switch (tribology.filmThickness) {
	case BS_In::Tribology::FilmThickness::FilmThicknessNothing:
		this->TB.FT = new Tribology::FilmThicknessNothing();
		break;
	case BS_In::Tribology::FilmThickness::HamrockDowsonHc:
		this->TB.FT = new Tribology::HamrockDowsonHc();
		break;
	case BS_In::Tribology::FilmThickness::HamrockDowsonHmin:
		this->TB.FT = new Tribology::HamrockDowsonHmin();
		break;
	}
	this->TB.TR = new Tribology::AiharaT();
	switch (tribology.hysteresis) {
	case BS_In::Tribology::Kakuta:
		this->TB.HY = new Tribology::Kakuta();
		this->TB.fh = tribology.hysteresis_factor;
		break;
	case BS_In::Tribology::HysteresisNothing:
		this->TB.HY = new Tribology::HysteresisNothing();
		this->TB.fh = 0;
	}

	return;
}

void BS_BallCylinderPair::init_Oil(const BS_In::Oil & oil) {

	this->LB.eta = oil.eta;		// �S�x�iPa*s�j
	this->LB.beta = oil.beta;	// ���x�S�x�W��
	this->LB.k = oil.k;			// ���M�`����
	this->LB.alpha = oil.alpha;	// ���͔S�x�W��
	this->LB.lm = oil.lm;		// ���j�X�J�X����

	return;
}



// �v�Z���ʂ������o�ϐ��Ɋi�[���郁�\�b�h�D
void BS_BallCylinderPair::save(
	int ig,
	const Vector3d&eta,
	const Vector3d&Fn,
	const Vector3d&Fs,
	const Vector3d&us,
	const Vector3d&ur,
	double dx,
	double cos_alp,
	double sin_alp,
	double a,
	double b,
	double h,
	const Vector3d&p,
	double lambda,
	double fratio,
	double Pmax
) {
	this->SV.eta = eta;
	this->SV.GV[ig].Fn = Fn;
	this->SV.GV[ig].Fs = Fs;
	this->SV.GV[ig].us = us;
	this->SV.GV[ig].ur = ur;
	this->SV.GV[ig].dx = dx;
	this->SV.GV[ig].phi = atan2(sin_alp, cos_alp);
	this->SV.GV[ig].a = a;
	this->SV.GV[ig].b = b;
	this->SV.GV[ig].h = h;
	this->SV.GV[ig].p = p;
	this->SV.GV[ig].lambda = lambda;
	this->SV.GV[ig].fratio = fratio;
	this->SV.GV[ig].Pmax = Pmax;
	return;
}

// ���݂̃{�[���ʒu�ɂ�����a�i�s�������擾���郁�\�b�h�i�������W�n�j
Vector3d BS_BallCylinderPair::get_e(void) {

	Vector3d BL_x = this->BL->x;
	Vector3d eta = this->CY->to_etacoord(iSP, BL_x);
	double th = eta[0];

	Matrix3d xyz2eta = this->CY->get_xyz2eta(iSP, th);
	Vector3d e = Vector3d(1.0, 0.0, 0.0);
	Vector3d e_ = this->CY->to_inertialvector(e, xyz2eta);

	return e_;
}

void BS_BallCylinderPair::save(BS_Out::BallCylinderPair & OUT) {

	Matrix3d xyz2eta = this->CY->get_xyz2eta(iSP, this->SV.eta[0]);

	for (int k = 0; k < 3; k++)
		OUT.eta[k] = this->SV.eta[k];
	for (int j = 0; j < 2; j++) {
		Vector3d Fn_ = this->CY->to_etavector(this->SV.GV[j].Fn, xyz2eta);
		Vector3d Fs_ = this->CY->to_etavector(this->SV.GV[j].Fs, xyz2eta);
		Vector3d us_ = this->CY->to_etavector(this->SV.GV[j].us, xyz2eta);
		Vector3d ur_ = this->CY->to_etavector(this->SV.GV[j].ur, xyz2eta);
		for (int k = 0; k < 3; k++) {
			OUT.GV[j].Fn[k] = Fn_[k];
			OUT.GV[j].Fs[k] = Fs_[k];
			OUT.GV[j].us[k] = us_[k];
			OUT.GV[j].ur[k] = ur_[k];
		}
		OUT.GV[j].dx = this->SV.GV[j].dx;
		OUT.GV[j].phi = this->SV.GV[j].phi;
		OUT.GV[j].a = this->SV.GV[j].a;
		OUT.GV[j].b = this->SV.GV[j].b;
		OUT.GV[j].h = this->SV.GV[j].h;
		OUT.GV[j].fratio = this->SV.GV[j].fratio;
		OUT.GV[j].lambda = this->SV.GV[j].lambda;
		OUT.GV[j].Pmax = this->SV.GV[j].Pmax;
		write_slice(j, xyz2eta, OUT);
	}
	return;
}

// �X���C�X�Ђ̌��ʂ��o�͗p�\���̂Ɋi�[
void BS_BallCylinderPair::write_slice(
	int ig,						// in:	[-]: �a�ԍ�
	Matrix3d xyz2eta,			// in:	[-]: �a�������W�n�ϊ��s��
	BS_Out::BallCylinderPair&OUT	// out:    : ��-���O�֐ڐG�v�Z�o��
) {
	for (int j = 0; j < this->GV[ig].n; j++) {
		OUT.GV[ig].SL[j].f_arr = this->SV.GV[ig].SL[j].f_arr;
		OUT.GV[ig].SL[j].mu_cl = this->SV.GV[ig].SL[j].mu_cl;
		OUT.GV[ig].SL[j].mu_tr = this->SV.GV[ig].SL[j].mu_tr;
		Vector3d fs_ = this->CY->to_etavector(this->SV.GV[ig].SL[j].fs, xyz2eta);
		Vector3d ts_ = this->CY->to_etavector(this->SV.GV[ig].SL[j].ts, xyz2eta);
		Vector3d us_ = this->CY->to_etavector(this->SV.GV[ig].SL[j].us, xyz2eta);
		Vector3d p_ = this->CY->to_etavector(this->SV.GV[ig].SL[j].ps, xyz2eta);

		for (int k = 0; k < 3; k++) {
			OUT.GV[ig].SL[j].fs_[k] = fs_[k];
			OUT.GV[ig].SL[j].ts_[k] = ts_[k];
			OUT.GV[ig].SL[j].us_[k] = us_[k];
			OUT.GV[ig].SL[j].ps_[k] = p_[k];
		}
	}
	return;
}


BS_BallCylinderPair::BS_BallCylinderPair() {
	for (int i = 0; i < 2; i++) {
		this->GV[i].ps = NULL;
		this->GV[i].r = NULL;
	}
	return;
}

BS_BallCylinderPair::~BS_BallCylinderPair() {
	for (int i = 0; i < 2; i++) {
		if (this->GV[i].ps != NULL)
			delete[] this->GV[i].ps;
		if (this->GV[i].r != NULL)
			delete[] this->GV[i].r;
	}
	return;
}


//int is = this->is;
//double th = this->eta[0];
//Matrix3d xyz2eta = this->CY->get_xyz2eta(is, th);
//Vector3d xai = this->CY->to_inertialvector(Vector3d::UnitX(), xyz2eta);

//for (int ig = 0; ig < 2; ig++) {
//	Vector3d p  = this->GV[ig].p;
//	double Fn = this->GV[ig].F;
//	double a  = this->GV[ig].a;

//	int ms = this->GV[ig].n;
//	VectorXd F_arr = Tribology::FroceSlice(Fn, ms);

//	this->CY->calc_slice(is, ig, th, p, xai, a, ms, this->GV[ig].ps);
//	double ur = this->get_ur(p).norm();

//	for (int is = 0; is < ms; is++) {
//		Vector3d ps = this->GV[ig].ps[is];
//		Vector3d us = this->get_us(ps);
//		double us_norm =  us.norm();
//		double mu = Tribology::CoulombsModel(this->GV[ig].mu, us_norm, ur);
//		Vector3d F_ = -us * F_arr[is] / us_norm * mu;
//		Vector3d T_ = this->BL->calc_Torque(p, F_);

//		F += F_;
//		T += T_;
//	}
//}
//void BS_BallCylinderPair::get_vFvTsum0(Vector3d & vF, Vector3d & vT) {
//
//	vF = vT = Vector3d::Zero();
//
//	for (int ig = 0; ig < 2; ig++) {
//
//		Vector3d p = this->GV[ig].p;
//		Vector3d v = this->get_us(p);
//		double   F = this->GV[ig].F;
//		Vector3d vF_= v * F;
//		Vector3d vT_= this->BL->calc_Torque(p, vF_);
//		vF += vF_;
//		vT += vT_;
//	}
//	return;
//}
	//F = T = Vector3d::Zero();

	//int is = this->is;
	//double th = this->eta[0];
	//Matrix3d xyz2eta = this->CY->get_xyz2eta(is, th);
	//Vector3d xai = this->CY->to_inertialvector(Vector3d::UnitX(), xyz2eta);

	//for (int ig = 0; ig < 2; ig++) {
	//	Vector3d p  = this->GV[ig].p;
	//	double Fn = this->GV[ig].F;
	//	double a  = this->GV[ig].a;

	//	int ms = this->GV[ig].n;
	//	Tribology::ForceSlice(Fn, ms, this->GV[ig].fs);

	//	this->CY->calc_slice(is, ig, th, p, xai, a, ms, this->GV[ig].ps);
	//	double ur = this->get_ur(p).norm();

	//	for (int is = 0; is < ms; is++) {
	//		Vector3d ps = this->GV[ig].ps[is];
	//		Vector3d us = this->get_us(ps);
	//		double us_norm =  us.norm();
	//		double mu = this->GV[ig].mu;	// Tribology::CoulombsModel(this->GV[ig].mu, us_norm, ur);
	//		Vector3d F_ = -us * this->GV[ig].fs[is] / us_norm * mu;
	//		Vector3d T_ = this->BL->calc_Torque(p, F_);

	//		F += F_;
	//		T += T_;
	//	}
	//}
//// �׏d�~�a�\�ʑ��x�̑��a���擾���郁�\�b�h�i�������W�n�j
//Vector3d BS_BallCylinderPair::get_vFsum0(void) {
//
//	Vector3d vF(0.0, 0.0, 0.0);
//
//	for (int ig = 0; ig < 2; ig++) {
//
//		Vector3d p = this->GV[ig].p;
//		Vector3d v = this->CY->surface_velocity(p);
//		double   F = this->GV[ig].F;
//
//		vF += v * F;
//	}
//	return vF;
//}