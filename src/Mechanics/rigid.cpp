#include "Rigid.h"

double Rigid::l;
double Rigid::t;
Vector3d Rigid::g;

// �R���X�g���N�^�C�������x�N�g���i�萔�l�j���`
Rigid::Rigid(void) {
	// ������������x�����������Ƃ���D
	this->ax0 << 1.0, 0.0, 0.0;
	this->ay0 << 0.0, 1.0, 0.0;
	this->az0 << 0.0, 0.0, 1.0;
	
	this->x_const = Eigen::Array<bool, 3, 1>(false, false, false);
	this->Rx_const = Eigen::Array<bool, 3, 1>(false, false, false);

	return;
}

// �l�̐ݒ胁�\�b�h�D
void Rigid::set_param(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w) {
	this->x = x;
	this->v = v;
	this->q = q;
	this->w = w;
	this->q.normalize();
	return;
}

// �l�̎擾���\�b�h�D
void Rigid::get_param(Vector3d&x, Vector3d&v, Quaterniond&q, Vector3d&w) {
	x = this->x;
	v = this->v;
	q = this->q;
	w = this->w;
	return;
}


// �l�̐ݒ胁�\�b�h�D�������ʂœn���ꂽ�l��L�����ɂ��ăp�����^�Ƃ��đ������D
void Rigid::set_y(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w) {
	this->x = x * Rigid::l;
	this->v = v * Rigid::l / Rigid::t;
	this->q = q;
	this->w = w / Rigid::t;
	this->q.normalize();
	return;
}


// �l�̎擾���\�b�h�D�L�����ʂŕێ����Ă���l�𖳎����ʂɂ��ĕԂ����\�b�h�D
void Rigid::get_y(Vector3d&x, Vector3d&v, Quaterniond&q, Vector3d&w) {
	x = this->x / Rigid::l;
	v = this->v / Rigid::l * Rigid::t;
	q = this->q;
	w = this->w * Rigid::t;
	return;
}

// ���̂܂�܂����ǌ��݂̑��x���o�͂���D
Vector3d Rigid::get_dxdt() {
	return this->v;
}

// ���ʂ���������x�ω��ʂ��o�͂���D
Vector3d Rigid::get_dvdt(const Vector3d&F) {
	Vector3d dvdt = this->m_inv * F;
	return dvdt;
}

// ���݂̃N�H�[�^�j�I������p���x�̓��͂ɂ������N�H�[�^�j�I���ω��ʂ��o�͂���D
// �N�H�[�^�j�I�����m�̊O�ςƂ��Čv�Z
Quaterniond Rigid::get_dqdt(void) {
	Quaterniond w_half(0.0, 0.5*this->w.x(), 0.5*this->w.y(), 0.5*this->w.z());
	Quaterniond dqdt = w_half * this->q;
	return dqdt;
}

// ���݂̊p���x����g���N�̓��͂ɂ������p���x�ω��ʂ��o�͂���D
Vector3d Rigid::get_dwdt(const Vector3d&T) {
	Vector3d Iw = this->I.cwiseProduct(this->w);
	Vector3d wIw = this->w.cross(Iw);
	Vector3d dwdt = this->I_inv.cwiseProduct(T - wIw);

	return dwdt;
}

// F, T �̓��͂ɂ�� x, v, q, w �S�Ă̎��Ԕ����l��Ԃ����\�b�h�D�S�������t���D�l��S�Ė������ʂɂ��ĕԂ��D
void Rigid::get_dydt(
	const Vector3d&F,	// in:	[N]:	�O���׏d�D
	const Vector3d&T,	// in:	[Nm]:	�O���g���N�D
	double*dydt			// out:	[-]:	�S�����l�i�������ʁj�D
) {
	Vector3d	dxdt = this->get_dxdt() / Rigid::l * Rigid::t;
	Quaterniond	dqdt = this->get_dqdt();
	Vector3d dvdt = this->get_dvdt(F) / Rigid::l * Rigid::t * Rigid::t; 
	Vector3d dwdt = this->get_dwdt(T) * Rigid::t * Rigid::t;

	dvdt = (this->x_const == true).select(Vector3d::Zero(), dvdt); // ���N�H�[�^�j�I���ɍS�������������Ă��Ȃ����߁C�������S������Ȃ�
	dwdt = (this->Rx_const == true).select(Vector3d::Zero(), dwdt);


	dydt[0]  = dxdt.x(); dydt[1]  = dxdt.y(); dydt[2]  = dxdt.z();
	dydt[3]  = dvdt.x(); dydt[4]  = dvdt.y(); dydt[5]  = dvdt.z();
	dydt[6]  = dqdt.w() * Rigid::t; dydt[7]  = dqdt.x() * Rigid::t; dydt[8]  = dqdt.y() * Rigid::t; dydt[9]  = dqdt.z() * Rigid::t;
	dydt[10] = dwdt.x(); dydt[11] = dwdt.y(); dydt[12] = dwdt.z();
	return;
}

// �S��������`
void Rigid::set_const(const bool (&v_const)[3], const bool (&w_const)[3]) {
	for (int i = 0; i < 3; i++) {
		this->x_const[i] = v_const[i];
		this->Rx_const[i] = w_const[i];
	}
	return;
}

// ���͂��ꂽ�����n���W���C�����̌n�̍��W�ɕϊ����郁�\�b�h�D
Vector3d Rigid::to_mycoord(const Vector3d&x) {
	Quaterniond qinv = this->q.inverse();
	Vector3d dx = x - this->x;
	Vector3d y  = qinv * dx;
	return y;
}

// ���͂��ꂽ�����̌n�̍��W���C�����n���W�ɕϊ����郁�\�b�h�D
Vector3d Rigid::to_inecoord(const Vector3d&y) {
	Vector3d dx = this->q * y;
	Vector3d x  = dx + this->x;
	return x;
}

// ���͂��ꂽ�����n���x���C�����̌n�̑��x�ɕϊ����郁�\�b�h�D
Vector3d Rigid::to_myvelocity(const Vector3d&v) {
	Quaterniond qinv = this->q.inverse();
	Vector3d dv = v - this->v;
	Vector3d V  = qinv * dv;
	return V;
}

// ���͂��ꂽ�����̌n�̍��W���C�����n���W�ɕϊ����郁�\�b�h�D
Vector3d Rigid::to_inevelocity(const Vector3d&V) {
	Vector3d dv = this->q * V;
	Vector3d v  = dv + this->v;
	return v;
}

// ���͂��ꂽ�����n�x�N�g���������̌n�̃x�N�g���ɕϊ����郁�\�b�h�D
Vector3d Rigid::to_myvector(const Vector3d&x) {
	Quaterniond qinv = this->q.inverse();
	Vector3d y  = qinv * x;
	return y;
}

// ���͂��ꂽ�����̌n�̃x�N�g���������n�x�N�g���ɕϊ����郁�\�b�h�D
Vector3d Rigid::to_inevector(const Vector3d&y) {
	Vector3d x  = this->q * y;
	return x;
}

// ���͂��ꂽ�������W�n��̓_�ɂ��āC���̑��x�����߂�i�����̑̓��̓_���͔��肵�Ȃ��j
Vector3d Rigid::surface_velocity(const Vector3d&x) {
	Vector3d r = x - this->x;
	Vector3d V = w.cross(r);
	V = V + this->v;
	return V;
}


// ���͂��ꂽ�������W�n��̓_�ɂ�����͂Ńg���N�����߂�D�i�����̑̓��̓_���͔��肵�Ȃ��j
Vector3d Rigid::calc_Torque(const Vector3d&x, const Vector3d&F) {
	Vector3d r = x - this->x;
	Vector3d T = r.cross(F);
	return T;
}


// ���͂��ꂽ�ʒu�E�����ɗ͂����������ۂ̃g���N�̌��������߂�D�i�����̑̓��̓_���͔��肵�Ȃ��j
Vector3d Rigid::calc_TorqueDirection(const Vector3d&x, const Vector3d&u) {
	Vector3d r = x - this->x;
	Vector3d T = r.cross(u);
	return T.normalized();
}


// ���݂̎��������Z�o���郁�\�b�h�D�i�������W�n�C�m����=1�j
Vector3d Rigid::get_ax(void) {
	Vector3d ax = this->to_inevector(this->ax0);
	return ax;
}


// ���݂̕��ʊp180���̌a�������Z�o���郁�\�b�h�D����D�i�������W�n�C�m����=1�j
Vector3d Rigid::get_ay(void) {
	Vector3d ay = this->to_inevector(this->ay0);
	return ay;
}


// ���݂̕��ʊp270���̌a�������Z�o���郁�\�b�h�D����D�i�������W�n�C�m����=1�j
Vector3d Rigid::get_az(void) {
	Vector3d az = this->to_inevector(this->az0);
	return az;
}


// ���͂��ꂽ�d�͉����x�x�N�g�����玩�g�̎��ʂ̐ς��Ƃ�C�׏d�x�N�g�����Z�o����D
Vector3d Rigid::get_mg(void) {
	return this->m * Rigid::g;
}

// ���͂��ꂽ�x�N�g���́C�������������������C���x�N�g���ƒ��p�ȃx�N�g���ɐ��`����D
Vector3d Rigid::remove_ax(const Vector3d&x) {
	Vector3d ax = this->get_ax();
	return x - ax.dot(x) * ax;
}

//�@���͂����������������悤�CQuaternion���v�Z���ݒ肷��Dq.x()��0�ŌŒ�D
void Rigid::set_ax(const Vector3d&ax) {
	Vector3d ax_ = ax.normalized();
	double y2z2 = 0.5 * (1.0 - ax_.x());
	double sq_1_y2z2 = sqrt(1 - y2z2);
	double y = -0.5 * ax_.z() / sq_1_y2z2;
	double z =  0.5 * ax_.y() / sq_1_y2z2;

	Quaterniond q_= Quaterniond(sq_1_y2z2, 0.0, y, z);
	this->q = q_;

}

// �g���N�Ɖ׏d�������o�ϐ��ɑ��
void Rigid::set_FT(const Vector3d & F, const Vector3d & T) {
	this->F = F;
	this->T = T;
	return;
}

void Rigid::set_mI(double m, const Vector3d & I) {
	this->m = m;
	this->m_inv = 1.0 / m;
	this->I = I;
	this->I_inv = I.cwiseInverse();
	return;
}

// x,v,q,w ��4�p�����[�^���O�ɏ����o�����\�b�h�D�i13�ϐ��j
void Rigid::save(double*x, double*v, double*q, double*w, double*ax, double*F, double*T) {

	Vector3d ax_ = this->get_ax();
	for (int i = 0; i < 3; i++) {
		x[i] = this->x[i];
		v[i] = this->v[i];
		w[i] = this->w[i];
		ax[i] = ax_[i];
		F[i] = this->F[i];
		T[i] = this->T[i];
	}
	q[0] = this->q.x();
	q[1] = this->q.y();
	q[2] = this->q.z();
	q[3] = this->q.w();
	return;
}

Rigid::~Rigid(void) {
}
