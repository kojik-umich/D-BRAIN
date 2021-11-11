#include "B4P_BallCagePair.h"

// �ʁE�ێ���I�u�W�F�N�g�������o�ϐ��ɑ��
void B4P_BallCagePair::link(Ball*BL, bal_Cage*CG, int np) {

	this->BL = BL;
	this->CG = CG;
	this->np = np;

	return;
}

// �p�����[�^�������o�ϐ��ɑ��
void B4P_BallCagePair::init(const B4P_In&FI) {

	this->mu = FI.BCP.mu;
	this->zeta = FI.BCP.dzeta;
	if (false) {
		this->DF = new Tribology::Tsuji();
		Tribology::BrewHamrock BH;
		double Rx =  1. / (1. / FI.Snap.R + 1. / this->BL->r);
		double E = 1. / (1. / FI.Cage.E + 1. / FI.BL[0].E);
		double k, a, b;
		BH.calc(Rx, Rx, 1, E, k, a, b);
		printf("%f\n", k);
	}
	else {
		this->DF = new Tribology::KelvinVoigt();
	}
	this->m = Tribology::ReducedMass(this->BL->m, this->CG->m);

	return;
}

// ��-�ێ���Ԃɓ����׏d�E�g���N���v�Z
void B4P_BallCagePair::calc_force(Vector3d&Fbc, Vector3d&Tbc, Vector3d&Fcb, Vector3d&Tcb) {

	// �߂�l�̐錾����я������D
	double k[_MAX_CONTACT_], dx[_MAX_CONTACT_]; 
	int ptt[_MAX_CONTACT_];	// �ڐG�p�^�[��
	Vector3d p[_MAX_CONTACT_], Fn[_MAX_CONTACT_], Fs[_MAX_CONTACT_];	// �ڐG�_�ʒu�ɂ����鍄���D
	for (int i = 0; i < _MAX_CONTACT_; i++) {
		p[i] = Fn[i] = Fs[i] = Vector3d::Zero();
		ptt[i] = 0;
	}
	// �ڐG�_�ʒu��ێ���ɔ��ʂ��Ă��炤�D
	int num_cont = this->CG->get_ContactPoint(this->BL->x, this->BL->r, this->np, p, k, ptt);

	// ���[�v�ŐώZ�����߂Ă������C���̑O�ɏ��������s���D
	Fbc = Tbc = Fcb = Tcb = Vector3d::Zero();

	// �ڐG�_�̐������׏d�E���C�v�Z�����C���Z���Ă����D
	for (int i = 0; i < num_cont; i++) {

		// �ڐG�׏d�̌v�Z
		double _dx = 0, dv = 0;
		Vector3d edir = Vector3d::Zero();
		Vector3d v = this->CG->surface_velocity(p[i]);
		bool c = this->BL->calc_Contact(p[i], v, _dx ,dv, edir);

		// �ڐG���Ȃ��ꍇ�͌�̏������ȗ��D
		if (!c) 
			continue;

		// �������l�����������׏d�̌v�Z
		double fn = this->DF->calc(k[i], this->zeta, this->m, dv, _dx);
		Vector3d _Fn = fn * -edir;

		// ���薀�C�̌v�Z
		Vector3d us = this->get_us(p[i]);
		double us_norm = us.norm();
		Vector3d _Fs = Vector3d::Zero();
		if (us_norm != 0.0) 
			_Fs = this->mu * _Fn.norm() * -us / us_norm;
		
		Fn[i] = _Fn;
		Fs[i] = _Fs;
		dx[i] = _dx;
		Fbc += _Fn + _Fs;
		Tbc += this->BL->calc_Torque(p[i], Fbc);
		Tcb += this->CG->calc_Torque(p[i], -Fbc);
	}
	// �v�Z���ʂ������o�ϐ��Ɋi�[
	for (int i = 0; i < _MAX_CONTACT_; i++) {
		this->p[i] = p[i];
		this->Fn[i] = Fn[i];
		this->Fs[i] = Fs[i];
		this->dx[i] = dx[i];
		this->ptt[i] = ptt[i];
	}
	Fcb = -Fbc;

	return;
}


// ���葬�x���Z�o���郁�\�b�h�i�������W�n�j
Vector3d B4P_BallCagePair::get_us(
	const Vector3d&p		// �ڐG�_�ʒu�D�i�������W�n�j
) {
	Vector3d ub = this->BL->surface_velocity(p);
	Vector3d uc = this->CG->surface_velocity(p);
	Vector3d us = ub - uc;
	Vector3d us_ = this->BL->remove_Normal(p, us);
	return us_;
}

// �ŐV�̏�Ԃ́C�ʂɂ�����C�ڐG�_�ʒu�C�ڐG�׏d�x�N�g���C���薀�C�̓x�N�g���C�]���薀�C�̓x�N�g���C���葬�x�x�N�g���C�]���葬�x�x�N�g���C�����_�l�C�׏d�x�������C�ڐG�p�̗]���E�����C�ڐG�ȉ~�����a�E�Z���a
void B4P_BallCagePair::save(B4P_Out::BallCagePair & BCP) {

	for (int i = 0; i < _MAX_CONTACT_; i++) {
		Vector3d p_ = this->CG->to_mycoord(this->p[i]);
		Vector3d Fn_ = this->CG->to_myvector(this->Fn[i]);
		Vector3d Fs_ = this->CG->to_myvector(this->Fs[i]);

		for (int j = 0; j < 3; j++) {
			BCP.CP[i].p[j] = this->p[i][j];
			BCP.CP[i].Fn[j] = this->Fn[i][j];
			BCP.CP[i].Fs[j] = this->Fs[i][j];
			BCP.CP[i].p_[j] = p_[j];
			BCP.CP[i].Fn_[j] = Fn_[j];
			BCP.CP[i].Fs_[j] = Fs_[j];
		}
		BCP.CP[i].dx = this->dx[i];
		BCP.CP[i].ptt = this->ptt[i];
	}
	return;
}

