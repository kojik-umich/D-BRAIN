/*******************************************************************************
!								"BS_Shaft.cpp"
!													2019/03/22	[Core-T]	���
!	�˂����i�ȉ��V���t�g�j�̓��f�����O�̓s����C
!	1. ���ޒP�̂̃V���t�g�i�p�[�c�j
!	2. �P�̃V���t�g�𕡐��g�ݍ��킹�����́i�A�Z���u���j
!	�̓��ނ���D�i����̃A�Z���u���V���t�g��1�̑g�ݍ��킹�����l���ł��Ȃ��j
!	���̃N���X�ł̓A�Z���u���V���t�g���`���Ă���C
!	Cylinder�N���X�ł˂��a�̏ڍׂȌ`����`���Ă���D
!
!*******************************************************************************/
#include "BS_Shaft.h"

void BS_Shaft::Nothing(void) {
}

void BS_Shaft::allocate(const std::vector<BS_In::Cylinder>&cylinders) {

	this->nCY = 1;
	this->CY = new BS_Cylinder[this->nCY];

	this->CY[0].allocate(cylinders[0].spiral.size());

	return;
}

void BS_Shaft::init(const std::vector<BS_In::Cylinder>&cylinders, const bool(&v_const)[3], const bool(&w_const)[3], const double(&ax0)[3], double v0, double w0) {

	this->CY[0].init(cylinders[0], v_const, w_const);

	this->rho = this->CY[0].rho = cylinders[0].density;
	this->nu = this->CY[0].nu = cylinders[0].poisson;
	this->E = this->CY[0].E = cylinders[0].young;

	this->w = this->get_ax() * w0;
	this->v = this->get_ax() * v0;

	Vector3d ax0_ = Vector3d(ax0) / ax0[0];		// x������1�ɂȂ�悤�Ɉꉞ�K�i���D
	this->tan_thy =-ax0_[2];
	this->tan_thz = ax0_[1];

	this->set_mI(cylinders[0].m, Vector3d(cylinders[0].Ix, cylinders[0].Iyz, cylinders[0].Iyz));

	this->set_const(v_const, w_const);

	//this->free_num = 6;
	//for (int i = 0; i < 3; i++) {
	//	this->free_num -= int(this->x_const[i]);
	//	this->free_num -= int(this->Rx_const[i]);
	//}
	return;
}

// �ʒu�̏������D�Ƃ肠�������_�ɒu���΂悢�H
void BS_Shaft::init_pos(double v0, double w0) {

	this->x = Vector3d::Zero();

	Vector3d ax = Vector3d::UnitX();
	this->set_ax(ax);
	this->v = v0 * ax;
	this->w = w0 * ax;

	this->set_dx();

	return;
}

// step0�i�����v�Z�j�p�̃C���^�t�F�C�X�D
void BS_Shaft::get_y0(double*y0) {

	y0[0] = this->x.x() / Rigid::l;
	y0[1] = this->x.y() / Rigid::l;
	y0[2] = this->x.z() / Rigid::l;

	y0[3] = this->q.y();
	y0[4] = this->q.z();

	return;
}

// step0�i�����v�Z�j�p�̃C���^�t�F�C�X�D(�S�������͗v�C��)
void BS_Shaft::set_y0(const double*y0, double v0, double w0) {

	Vector3d x_now = this->x;
	this->x = this->x_const.select(
		x_now,
		Vector3d(y0[0], y0[1], y0[2]) * Rigid::l
	);
	Vector3d q_ = this->Rx_const.select(
		Vector3d(this->q.x(), this->q.y(), this->q.z()),
		Vector3d(0.0, y0[3], y0[4])
	);
	this->q.y() = q_.y();
	this->q.z() = q_.z();
	this->q.normalize();

	Vector3d ax = this->get_ax();
	this->v = v0 * ax;
	this->w = w0 * ax;
	this->set_dx();

	return;
}

// �V���t�g�i�p�[�c�j�̐��l���V���t�g�i�A�Z���u���j�̐��l����ɍX�V�D�V���t�g�i�A�Z���u���j�ƃV���t�g�i�p�[�c�j�͓���ɓ������̂Ƃ���D
void BS_Shaft::set_dx(void) {
	this->CY[0].set_param(this->x, this->v, this->q, this->w);
	return;
}

// �l�̐ݒ胁�\�b�h�D�������ʂœn���ꂽ�l��L�����ɂ��ăp�����^�Ƃ��đ������D
void BS_Shaft::set_y_(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w) {

	this->set_y(x, v, q, w);
	this->set_dx();

	return;
}

// step0�i�����v�Z�j�p�̃C���^�t�F�C�X�D
void BS_Shaft::init_dyn0(void) {

	this->v = Vector3d::Zero();
	this->w = Vector3d::Zero();

	return;
}

void BS_Shaft::deinit_dyn0(double v0, double w0) {

	this->v = v0 * this->get_ax();
	this->w = w0 * this->get_ax();

	return;
}

// step0�i�����v�Z�j�p�̃C���^�t�F�C�X�D
void BS_Shaft::get_dyn_y0(double*y0) {

	for (int i = 0; i < 3; i++) {
		y0[i + 0] = this->x[i] / Rigid::l;
		y0[i + 3] = this->v[i] / Rigid::l * Rigid::t;
	}
	y0[6] = this->q.y() * 2;
	y0[7] = this->q.z() * 2;
	y0[8] = this->w.y() * Rigid::t;
	y0[9] = this->w.z() * Rigid::t;

	return;
}

// step0�i�����v�Z�j�p�̃C���^�t�F�C�X�D(�S�������͗v�C��)
void BS_Shaft::set_dyn_y0(const double*y0) {

	Vector3d x_now = this->x;
	Vector3d x = this->x_const.select(x_now,
		Vector3d(y0[0], y0[1], y0[2]) * Rigid::l);
	Vector3d v = Vector3d(y0[3], y0[4], y0[5])
		* Rigid::l / Rigid::t;
	Vector3d w = Vector3d(0, y0[8], y0[9])
		/ Rigid::t;

	Vector3d t_now = 2 * Vector3d(0, this->q.y(), this->q.z());
	Vector3d t = this->Rx_const.select(t_now, Vector3d(0, y0[6], y0[7]));
	double   qy = 0.5 * t.y();
	double   qz = 0.5 * t.z();
	double   qw = sqrt(1.0 - qy * qy - qz * qz);

	this->x = x;
	this->v = v;
	this->q = Quaterniond(qw, 0, qy, qz);
	this->w = w;

	this->set_dx();

	return;
}

// F, T �̓��͂ɂ�� x, v, q, w �S�Ă̎��Ԕ����l��Ԃ����\�b�h�D�l��S�Ė������ʂɂ��ĕԂ��D
void BS_Shaft::get_dydt_(
	const Vector3d&F,	// in:	[N]:		�O���׏d�D
	const Vector3d&T,	// in:	[Nm]:		�O���g���N�D
	double dvdt0,		// in:	[m/s^2]:	���x�̎��ԕω���
	double dwdt0,		// in:	[rad/s^2]	��]���̎��ԕω���
	double*dydt			// out:	[any]:		�S�����l�D
) {
	// �p���S���i�p���x�ƃN�H�[�^�j�I���̗������S���C������x���̉�]�̂ݕ��̍��W�n�ōS���j
	Vector4d dqdt = this->get_dqdt_(this->Rx_const.y(), this->Rx_const.z(), this->tan_thy, this->tan_thz) * Rigid::t;
	Vector3d dwdt = this->get_dwdt_(T, this->Rx_const.y(), this->Rx_const.z(), dwdt0) * Rigid::t * Rigid::t;

	// �ψʍS��(�������W�n�ōS��)
	Vector3d dxdt = this->get_dxdt() / Rigid::l * Rigid::t;
	Vector3d _dvdt = this->get_dvdt(F) / Rigid::l * Rigid::t * Rigid::t;
	Vector3d dvdt = (this->x_const == true).select(Vector3d::Zero(), _dvdt);
	if (this->x_const.x()) {
		dvdt += Vector3d(dvdt0, 0, 0);
	}

	dydt[0] = dxdt.x(); dydt[1] = dxdt.y(); dydt[2] = dxdt.z();
	dydt[3] = dvdt.x(); dydt[4] = dvdt.y(); dydt[5] = dvdt.z();
	dydt[6] = dqdt.w(); dydt[7] = dqdt.x(); dydt[8] = dqdt.y(); dydt[9] = dqdt.z();
	dydt[10] = dwdt.x(); dydt[11] = dwdt.y(); dydt[12] = dwdt.z();

	return;
}

// �p�x�S���t���N�H�[�^�j�I���̎��Ԕ������Dy/z���ɉ�]�S����������ꍇ�C�v�Z�덷�Ŏp�������X�ɕω����邽�߁C�S���������Ă���D
Vector4d BS_Shaft::get_dqdt_(bool wy_const, bool wz_const, double tan_thy, double tan_thz) {
	// Quaterniond �̃R���X�g���N�^��(w, x, y, z)�̏���
	Quaterniond w_half(0.0, 0.5*this->w.x(), 0.5*this->w.y(), 0.5*this->w.z());
	Vector4d qw_half = (w_half * q).coeffs();	// �N�H�[�^�j�I���͉��Z����`����Ă��Ȃ����߁C�x�N�g���z��Ōv�Z�D
	Vector4d dqdt;
	double tau = 0.1;			// �ɘa�W���i�{���͓��̓t�@�C��������͂ł���悤�ɂ���ׂ��j
	double qx = q.x(), qy = q.y(), qz = q.z(), qw = q.w();

	// ���S�S���͊p���x�S�������Ŗ��Ȃ����߁C�p���͍S���Ȃ�
	if (wy_const && wz_const) {
		dqdt = qw_half;
	}
	// wy�������S������Ă���ꍇ�Cwy�����Ǝ������ɉ����x�����������Ȃ��悤�ɂ���D
	// (= wy�P�ʃx�N�g���Ǝ������x�N�g���̊O�ς̌����݂̂ɉ�������悤�ɂ���)
	else if (wy_const && !wz_const) {
		// Vector4d �̃R���X�g���N�^��(x, y, z, w)�̏���
		Vector4d _Cq(qz + tan_thy * qx, -qw - tan_thy * qy,
			qx - tan_thy * qz, -qy + tan_thy * qw);
		Vector4d Cq = _Cq * 2;
		double C = 2 * (qx * qz - qy * qw) +
			(qx * qx - qy * qy - qz * qz + qw * qw) * tan_thy;

		// 0����̔���(��������Ȃ�)
		double Cq_square = Cq.squaredNorm();
		if (Cq_square < 1e-20) {
			dqdt = qw_half;
		}
		else {
			dqdt = qw_half - Cq * (qw_half.dot(Cq) + 1 / tau * C) / Cq_square;
		}
	}
	// wz�������S������Ă���ꍇ�Cwz�����Ǝ������ɉ����x�����������Ȃ��悤�ɂ���D
	// (= wy�P�ʃx�N�g���Ǝ������x�N�g���̊O�ς̌����݂̂ɉ�������悤�ɂ���)
	else if (!wy_const && wz_const) {
		// Vector4d �̃R���X�g���N�^��(x, y, z, w)�̏���
		Vector4d _Cq(qy - tan_thz * qx, qx + tan_thz * qy,
			qw + tan_thz * qz, qy - tan_thz * qw);
		Vector4d Cq = _Cq * 2;
		double C = 2 * (qx * qy + qz * qw) +
			(qx * qx - qy * qy - qz * qz + qw * qw) * tan_thz;

		// 0����̔���i��������Ȃ��j
		double Cq_sq = Cq.squaredNorm();
		if (Cq_sq < 1e-20) {
			dqdt = qw_half;
		}
		else {
			dqdt = qw_half - Cq * (qw_half.dot(Cq) + 1 / tau * C) / Cq_sq;
		}
	}
	// �������̉�]�̂ݍS���D
	else {
		dqdt = qw_half;
	}
	return dqdt;
}


// �p���x�̎��Ԕ������D�����]�����ƍS������Ă��鐬���͕ω����Ȃ��悤�ɐݒ�D
Vector3d BS_Shaft::get_dwdt_(const Vector3d&T, bool wy_const, bool wz_const, double dwdt0) {
	Vector3d dwdt;

	// ���S�S��
	if (wy_const && wz_const) {
		Vector3d ax = this->get_ax();
		dwdt = ax * dwdt0;
	}
	// wy�������S������Ă���ꍇ�Cwy�����Ǝ������ɉ����x�����������Ȃ��悤�ɂ���D
	// (= wy�P�ʃx�N�g���Ǝ������x�N�g���̊O�ς̌����݂̂ɉ�������悤�ɂ���)
	else if (wy_const && !wz_const) {
		Vector3d ax = this->get_ax();
		Vector3d ey = Vector3d(0, 1, 0);
		Vector3d ew = ax.cross(ey);			// y�����Ǝ������ɐ����ȃx�N�g��
		Vector3d _dwdt = this->get_dwdt(T);
		dwdt = _dwdt.dot(ew) * ew + ax * dwdt0;
	}
	// wz�������S������Ă���ꍇ�Cwz�����Ǝ������ɉ����x�����������Ȃ��悤�ɂ���D
	// (= wy�P�ʃx�N�g���Ǝ������x�N�g���̊O�ς̌����݂̂ɉ�������悤�ɂ���)
	else if (!wy_const && wz_const) {
		Vector3d ax = this->get_ax();
		Vector3d ez = Vector3d(0, 0, 1);
		Vector3d ew = ax.cross(ez);			// z�����Ǝ������ɐ����ȃx�N�g��
		Vector3d _dwdt = this->get_dwdt(T);
		dwdt = _dwdt.dot(ew) * ew + ax * dwdt0;
	}
	// �������̉�]�̂ݍS���D
	else {
		Vector3d ax = this->get_ax();
		Vector3d dwdt_ = this->get_dwdt(T);
		dwdt = this->remove_ax(dwdt_) + ax * dwdt0;
	}
	return dwdt;
}


void BS_Shaft::set_Fall_Zero(void) {

	this->F = this->T = Vector3d::Zero();

	for (int i = 0; i < this->nCY; i++)
		this->CY[i].F = this->CY[i].T = this->CY[i].Fs = this->CY[i].Ts = Vector3d::Zero();

	return;
}

BS_Shaft::BS_Shaft() {
	this->CY = NULL;
	return;
}

BS_Shaft::~BS_Shaft() {
	if (this->CY != NULL)
		delete[] this->CY;
	return;
}





