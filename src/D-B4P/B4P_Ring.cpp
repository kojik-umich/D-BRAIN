/******************************************************************************* 
!								"B4P_Ring.cpp"	
!													2019/03/22	[Core-T]	���
!	
!	
!*******************************************************************************/
#include "B4P_Ring.h"
#include<iostream>
void B4P_Ring::Nothing(void) {
}

// �ڐG�ȉ~���g�[���X�ɉ����ăX���C�X���郁�\�b�h�D(�֐����̍��W�ƃx�N�g���͂��ׂĊ������W�n)
void B4P_Ring::calc_slice
	(
	const Vector3d&p,	// in:	[m]:	�ڐG�_�ʒu�̒��S�D�i�������W�n�j
	double a,			// in:	[m]:	�ڐG�ȉ~���a�D
	int msmax,			// in:	[-]:	�ڐG�ȉ~�̃X���C�X���D
	int ig,				// in:	[-]:	�a�ԍ��D0 or 1;
	double bl_r,		// in:	[m]:	�ʔ��a
	Vector3d*ps		// out:	[m]:	�X���C�X���ꂽ�e�ڐG�_�ʒu�D�z��T�C�Y=msmax�D�i�������W�n�j
	)
{
	// ���O�֒��S����ڐG�_�֌������x�N�g�����擾
	Vector3d xp = p - this->x;

	// ���݂̎������̒P�ʃx�N�g�����Z�o
	Vector3d ax = this->get_ax();
	
	// �ڐG�_�ɂ�����������x�N�g���i�������ƃx�N�g��xp�ƒ�������P�ʃx�N�g���j
	Vector3d cs  = (xp.cross(ax)).normalized();
	
	// �ڐG�_�ɂ�����a�����x�N�g���i�������ɐ����j
	Vector3d er = ax.cross(cs);
	
	// �������W�n�ł́C�_p�����~�ʂ̒��SO�D�i=�aR���S�Ɠ���Ɖ���j
	Vector3d O = this->GV[ig].Rx * ax + this->GV[ig].Rr * er + this->x;

	// �ȉ~���X���C�X�����ۂ́C�����̊Ԋu�D
	double da = 2 * a / msmax;

	// �~�̌v�Z�̓R�X�g�������邽�߁C�ڐG�ȉ~���\�����������Ƃ���2���ߎ����s���D
	Vector3d op = p - O;
	double op_norm = op.norm();

	if (op_norm < GV[ig].r) {
		double _t = 0;
		std::cout << "err!!!" << op_norm << " " << GV[ig].r << std::endl;
	}
	Vector3d cn = op / op_norm; // �a�f�ʂɐ����ȃx�N�g��
	Vector3d ct = cn.cross(cs);	// �ڐG�f�ʂƕ��s����cs�Ɛ����ȕ����x�N�g��
	double Rm = 2 * bl_r * this->GV[ig].r / (bl_r + this->GV[ig].r);
	// �e�X���C�X�̒��S�_ pi �̍��W���Z�o
	for (int i = 0; i < msmax; i++) {
		double xi_ = i * da - a + 0.5 * da;		// p����Ƃ����_p_i��ct�x�N�g�������̋���
		double yi_ = - 0.5 * xi_ * xi_ / Rm;	// p����Ƃ����_p_i��cn�x�N�g�������̋���
		ps[i] = p + cn * yi_ + ct * xi_;
	}
	return;
}

// �p�����[�^�������o�ϐ��ɑ��
void B4P_Ring::init(const B4P_In::Ring&IN, double pcd) {
	this->pcd  = pcd;					// ��PCD
	for (int i = 0; i < 2; i++) {
		this->GV[i].Rx = IN.Rox[i];			// pcd ���猩���������ψ�(-x��)
		this->GV[i].r = IN.R[i];				// �a�a(-x��)
		this->GV[i].r_inv = 1.0 / IN.R[i];	// �a�a�̋ȗ�(-x��)
		this->GV[i].Rr = IN.Rod[i] / 2;		// ���O�֒��S����a���S�̋��� [m]
		this->GV[i].h = IN.hedge[i];			// �a������[m](�V�K�ǉ�)
	}
	this->sigmap   =IN.rms;					// �e��rms(�V�K�ǉ�)
	this->m			= IN.m;					// �d��
	this->m_inv		= 1.0 / this->m;
	this->E  = IN.E;							// ����
	this->nu = IN.por;
	this->I			= Vector3d(IN.Ix, IN.Iyz, IN.Iyz);
	this->I_inv		= Vector3d(1/IN.Ix, 1/IN.Iyz, 1/IN.Iyz);
	return;
}

// �S�������p�̕��̍��W�n�ɕϊ����邽�߂̉�]�s������߂�D
// �S�������p���W�n�Fx'���������������Ă�����W�n�D�������Cx'�����ɂ͉�]���Ȃ��D
Matrix3d B4P_Ring::get_Rth() {
	Vector3d ex_d = this->get_ax();		// ���̍��W�nx'�����i�������j
	Vector3d ex = Vector3d(1, 0, 0);	// �������W�nx����
	// �ړ��O�ƈړ���̊p�x�����
	double angle = acos(ex.dot(ex_d));
	//std::cout << ex_d << std::endl;
	// ��]��(�O�ς����߂Đ��K��)
	Vector3d axis = ex.cross(ex_d).normalized();
	// ��]�ʂƉ�]�������]�s��𐶐�
	Matrix3d Rth = Eigen::AngleAxisd(angle, axis).matrix();

	return Rth;
}



// F, T �̓��͂ɂ�� x, v, q, w �S�Ă̎��Ԕ����l��Ԃ����\�b�h�D�S�������t���D�l��S�Ė������ʂɂ��ĕԂ��D
// �N�H�[�^�j�I���̍S�������� Baumgarte �̕��@��K�p
void B4P_Ring::get_dydt_(
	const Vector3d&F,	// in:	[N]:	�O���׏d�D
	const Vector3d&T,	// in:	[Nm]:	�O���g���N�D
	double*dydt			// out:	[any]:	�S�����l�D
) {

	// �p���S���i�p���x�ƃN�H�[�^�j�I���̗����ɍS����������j
	Vector4d dqdt = this->get_dqdt_(this->Rx_const.y(), this->Rx_const.z()) * Rigid::t;
	Vector3d dwdt = this->get_dwdt_(T, this->Rx_const.y(), this->Rx_const.z()) * Rigid::t * Rigid::t;

	// �ψʍS��
	Vector3d dxdt = this->get_dxdt() / Rigid::l * Rigid::t;
	Vector3d dvdt = this->get_dvdt(F) / Rigid::l * Rigid::t * Rigid::t;
	dvdt = (this->x_const == true).select(Vector3d::Zero(), dvdt);

	dydt[0] = dxdt.x(); dydt[1] = dxdt.y(); dydt[2] = dxdt.z();
	dydt[3] = dvdt.x(); dydt[4] = dvdt.y(); dydt[5] = dvdt.z();
	dydt[6] = dqdt.w(); dydt[7] = dqdt.x(); dydt[8] = dqdt.y(); dydt[9] = dqdt.z();
	dydt[10] = dwdt.x(); dydt[11] = dwdt.y(); dydt[12] = dwdt.z();
	return;
}

// �p�x�S���t���N�H�[�^�j�I���̎��Ԕ������D�N�H�[�^�j�I���͉��Z����`����Ă��Ȃ����߁C�x�N�g���z��Ōv�Z�D
Vector4d B4P_Ring::get_dqdt_(bool wy_const, bool wz_const){
	// Quaterniond �̃R���X�g���N�^��(w, x, y, z)�̏���
	Quaterniond w_half(0.0, 0.5*this->w.x(), 0.5*this->w.y(), 0.5*this->w.z());
	Vector4d qw_half = (w_half * q).coeffs();
	Vector4d dqdt;
	double tan_beta = 0;
	double tan_gamma = 0;
	double tau = 0.1;
	

	// ���S�S��
	if (wy_const && wz_const) {
		dqdt = qw_half;
	}
	// wy�������S������Ă���ꍇ�Cwy�����Ǝ������ɉ����x�����������Ȃ��悤�ɂ���D
	// (= wy�P�ʃx�N�g���Ǝ������x�N�g���̊O�ς̌����݂̂ɉ�������悤�ɂ���)
	else if (wy_const && !wz_const) {
		// Vector4d �̃R���X�g���N�^��(x, y, z, w)�̏���
		Vector4d _Cq(q.z() + tan_beta * q.x(), -q.w() - tan_beta * q.y(),
			q.x() - tan_beta * q.z(), -q.y() + tan_beta * q.w());
		Vector4d Cq = _Cq * 2;
		double C = 2 * (q.x() * q.z() - q.y() * q.w()) +
			(q.x() * q.x() - q.y() * q.y() - q.z() * q.z() + q.w() * q.w()) * tan_beta;
		// 0����̔���
		if (Cq.norm() < 1e-20) {
			dqdt = qw_half;
		}
		else {
			dqdt = qw_half - Cq * (qw_half.dot(Cq) + 1 / tau * C) / Cq.dot(Cq);
		}
	}
	// wz�������S������Ă���ꍇ�Cwz�����Ǝ������ɉ����x�����������Ȃ��悤�ɂ���D
	// (= wy�P�ʃx�N�g���Ǝ������x�N�g���̊O�ς̌����݂̂ɉ�������悤�ɂ���)
	else if (!wy_const && wz_const) {
		// Vector4d �̃R���X�g���N�^��(x, y, z, w)�̏���
		Vector4d _Cq(q.y() - tan_gamma * q.x(), q.x() + tan_gamma * q.y(),
			q.w() + tan_gamma * q.z(), q.y() - tan_gamma * q.w());
		Vector4d Cq = _Cq * 2;
		double C = 2 * (q.x() * q.y() + q.z() * q.w()) +
			(q.x() * q.x() - q.y() * q.y() - q.z() * q.z() + q.w() * q.w()) * tan_gamma;
		// 0����̔���
		if (Cq.norm() < 1e-20) {
			dqdt = qw_half;
		}
		else {
			dqdt = qw_half - Cq * (qw_half.dot(Cq) + 1 / tau * C) / Cq.dot(Cq);
		}
	}
	// �������̉�]�̂ݍS���D
	else {
		dqdt = qw_half;
	}
	return dqdt;
}

// �p���x�̎��Ԕ������D�����]�����ƍS������Ă��鐬���͕ω����Ȃ��悤�ɐݒ�D
Vector3d B4P_Ring::get_dwdt_(const Vector3d&T, bool wy_const, bool wz_const) {
	Vector3d dwdt;

	// ���S�S��
	if (wy_const && wz_const) {
		dwdt = Vector3d::Zero();
	}
	// wy�������S������Ă���ꍇ�Cwy�����Ǝ������ɉ����x�����������Ȃ��悤�ɂ���D
	// (= wy�P�ʃx�N�g���Ǝ������x�N�g���̊O�ς̌����݂̂ɉ�������悤�ɂ���)
	else if (wy_const && !wz_const) {
		Vector3d ax = this->get_ax();
		Vector3d ey = Vector3d(0, 1, 0);
		Vector3d ew = ax.cross(ey);
		Vector3d _dwdt = this->get_dwdt(T);
		dwdt = _dwdt.dot(ew) * ew;
	}
	// wz�������S������Ă���ꍇ�Cwz�����Ǝ������ɉ����x�����������Ȃ��悤�ɂ���D
	// (= wy�P�ʃx�N�g���Ǝ������x�N�g���̊O�ς̌����݂̂ɉ�������悤�ɂ���)
	else if (!wy_const && wz_const) {
		Vector3d ax = this->get_ax();
		Vector3d ez = Vector3d(0, 0, 1);
		Vector3d ew = ax.cross(ez);
		Vector3d _dwdt = this->get_dwdt(T);
		dwdt = _dwdt.dot(ew) * ew;
	}
	// �������̉�]�̂ݍS���D
	else {
		Vector3d dwdt_ = this->get_dwdt(T);
		dwdt = this->remove_ax(dwdt_);
	}
	return dwdt;
}

// �É�́@�����v�Z�F���ւ̍��W�E�p����z�� Xstf �Ɋi�[
void B4P_InnerRing::get_Xstf(double* Xstf) {

	Xstf[0] = this->x.x() / Rigid::l;
	Xstf[1] = this->x.y() / Rigid::l;
	Xstf[2] = this->x.z() / Rigid::l;

	Vector3d ax = this->get_ax();
	Xstf[3] = ax.y();
	Xstf[4] = ax.z();

	return;
}

// �É�́@�����v�Z�F���ւ̍��W�C�p�������
void B4P_InnerRing::set_Xstf(double* Xstf) {

	double x0 = sqrt(1.0 - Numeric::Square(Xstf[3]) - Numeric::Square(Xstf[4]));
	Vector3d ax(x0, Xstf[3], Xstf[4]);
	this->set_ax(ax);
	Vector3d ax_ = this->get_ax();

	this->x = Vector3d(Xstf[0], Xstf[1], Xstf[2]) * Rigid::l;

	return;
}

// �ʑ��p���犵�����W�n�̃x�N�g�����a���s�f�ʍ��W�n�֕ϊ����郁�\�b�h�D
Vector3d B4P_Ring::ine_to_XZvector
(							// out:	[any]:	0�����F�����O�ʑ��p�D1�����FX���W�D2�����FZ���W�D
	const Vector3d&a,		// in : [any]:	�ϊ��������x�N�g���D�i�������W�n�j
	const Matrix3d&xyz2XYZ	// in : [rad]:	�ʑ��p���̉�]�s��i�����O���W�n�j
) {
	Vector3d xyz = this->to_myvector(a);
	Vector3d XYZ = xyz2XYZ * xyz;
	return XYZ;
}

// �����O���W�n����a���s���W�n�֕ϊ�����s����擾���郁�\�b�h�D
Matrix3d B4P_Ring::get_xyz2XYZ
(						// out:	[any]:	0�����F�����O�ʑ��p�D1�����FX���W�D2�����FZ���W�D
	double th			// in : [rad]:	�ʑ��p�i�����O���W�n�j
) {
	double cos_th = cos(th);
	double sin_th = sin(th);
	Matrix3d xyz2XYZ;
	xyz2XYZ <<
		1.0, 0.0, 0.0,
		0.0, cos_th, sin_th,
		0.0, -sin_th, cos_th;
	return xyz2XYZ;
}

// �������W�n���烊���O�̍a�������W�n�ɕϊ����郁�\�b�h�D
Vector3d B4P_Ring::ine_to_XZcoord
(						// out:	[rad],[m]:	0�����F�����O�ʑ��p�D1�����FX���W�D2�����FZ���W�D
	const Vector3d& x	// in : [m]:		�������W�n�̕ϊ����������W�D
)
{
	Vector3d xyz = this->to_mycoord(x);
	double th = atan2(-xyz[1], xyz[2]);
	double X = xyz[0];
	double Z = sqrt(xyz[1] * xyz[1] + xyz[2] * xyz[2]) - 0.5 * this->pcd;
	Vector3d thXZ = Vector3d(th, X, Z);
	return thXZ;
}

// �����O�̍a�������W�n���犵�����W�n�ɕϊ����郁�\�b�h�D
Vector3d B4P_Ring::XZ_to_inecoord
(						// out:	[m]:		�������W�n�̍��W�D
	const Vector3d& thXZ	// in :	[rad],[m]:	�ϊ������������O�̍a�������W�n�D
)
{
	double th = thXZ[0];
	double x = thXZ[1];
	double R = 0.5 * this->pcd + thXZ[2];
	double y = -R * sin(th);
	double z = R * cos(th);
	Vector3d xyz = Vector3d(x, y, z);
	Vector3d xyz_ = this->to_inecoord(xyz);
	return xyz_;
}

// �^����ꂽ�C�ӂ̃x�N�g�����C����ʑ��p�ł̍a���p�f�ʏ��2�����x�N�g���Ɏˉe����D
Vector2d B4P_Ring::to_cross_sction
(					// out:	[Any]:	���͂Ɠ��������D��0�����F�������C��1�����F�a�����D
	const Vector3d&v,	// in:	[Any]:	�C�ӂ̃x�N�g���i�������W�n�j�D
	double th			// in:	[rad]:	�����ʒu�Ō����Ƃ��̃����O�ʑ��p�D
)
{
	Vector3d v_ = this->to_myvector(v);
	return Vector2d(v_[0], -sin(th)*v_[1] + cos(th)*v_[2]);
}


// ������������̉����畨�̂̊e�p�����[�^�����D
// �N�H�[�^�j�I���͐��l�덷�����邽�߁C���l������ɂ��S��������K�p�D
void B4P_Ring::set_y_(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w) {
	
	this->x = x * Rigid::l;
	this->v = v * Rigid::l / Rigid::t;
	this->w = w / Rigid::t;
	Quaterniond _q = q;
	// �N�H�[�^�j�I���͐��l�̌v�Z�덷�����邽�߁C�������x�N�g�����S������Ă�������ɓ����Ȃ��悤�ɂ���

	if (this->Rx_const[1] == true && _q.w() != 0) {
		_q.y() = _q.x() * _q.z() / _q.w();
	}
	if (this->Rx_const[2] == true && _q.w() != 0) {
		_q.z() = -_q.x() * _q.y() / _q.w();
	}
	this->q = _q;

	this->q.normalize();

}









//
//Vector3d B4P_Ring::to_myXYZ
//	(					// out:	[Any]:	���͂Ɠ��������D��0�����F�������C��1�����F�a�����D
//	const Vector3d&a	// in:	[Any]:	�C�ӂ̃x�N�g���i�������W�n�j�D
//	)
//{
//	Vector3d u_  = Vector3d(0.0, -sin(th), cos(th));	// �O�֍��W�n�ł̐i�s�����x�N�g���D
//	Vector3d u__ = this->to_inevelocity(u_);			// �������W�n�ւ̕ϊ��D
//	return u__;
//}
//
//void B4P_Ring::init(const B4P_In&FI, const InnerOuter&innerouter){
//	switch (innerouter)
//	{
//	case B4P_Ring::Inner:
//		this->pcd  = (FI.igrv_r_trs[0] + FI.igrv_r_trs[1]) / 2 * 2;
//		this->GV[0].x = -FI.X_disti/2;
//		this->GV[1].x = FI.X_disti/2;
//		this->GV[0].r = FI.rgi[0];
//		this->GV[1].r = FI.rgi[1];
//		this->GV[0].R = FI.igrv_r_trs[0];
//		this->GV[1].R = FI.igrv_r_trs[1];
//		this->m  = Numeric::pi/6*pcd*pcd*pcd*FI.rhop[1];
//		this->E  = FI.Ep[1];
//		this->nu = FI.etap[1];
//		this->m  = Numeric::pi/6*pcd*pcd*pcd*FI.rhop[1];	// �����O�����ǋ��̂ŉ���D
//
//	case B4P_Ring::Outer:
//		this->pcd  = (FI.ogrv_r_trs[0] + FI.ogrv_r_trs[1]) / 2 * 2;
//		this->GV[0].x = -FI.X_disto/2;
//		this->GV[1].x = FI.X_disto/2;
//		this->GV[0].r = FI.rgo[0];
//		this->GV[1].r = FI.rgo[1];
//		this->GV[0].R = FI.ogrv_r_trs[0];
//		this->GV[1].R = FI.ogrv_r_trs[1];
//		this->m  = Numeric::pi/6*pcd*pcd*pcd*FI.rhop[2];
//		this->E  = FI.Ep[2];
//		this->nu = FI.etap[2];
//		this->m  = Numeric::pi/6*pcd*pcd*pcd*FI.rhop[2];	// �����O�����ǋ��̂ŉ���D
//
//	default:
//		break;
//	}
//
//	double I_= this->m*pcd*pcd/10;
//	this->I << I_, I_, I_;
//	this->m_inv = 1.0 / this->m;
//	this->I_inv << 1 / I_, 1 / I_, 1 / I_;
//}
//$$Outer $$Inner[���O��],0���x,1�����O��,2por,3�e��rms,4�aR(-x),5�aR(+x),6���SO(-x),7���SO(+x),8�aR���S���a(-x),9�aR���S���a(+x),10������(-x), 11������(+x), 12��PCD
//void B4P_Ring::init(const VectorXd &IOringparam, double pcd){
//	//switch (innerouter)
//	//{
//	//case B4P_Ring::Inner:
//		//this->pcd  = IOringparam[12]; // ��PCD
//		this->pcd  = pcd; // ��PCD
//		this->GV[0].x = IOringparam[6];	// pcd ���猩���������ψ�(-x��)
//		this->GV[1].x = IOringparam[7];	// pcd ���猩���������ψ�(+x��)
//		this->GV[0].r = IOringparam[4];	// �a�a(-x��)
//		this->GV[1].r = IOringparam[5];	// �a�a(+x��)
//		this->GV[0].R = IOringparam[8];	// �a���Spcd [m]
//		this->GV[1].R = IOringparam[9];	// �a���Spcd [m]
//		this->GV[0].h = IOringparam[10];	// �a������[m](�V�K�ǉ�)
//		this->GV[1].h = IOringparam[11];	// �a������[m](�V�K�ǉ�)
//		this->sigmap   = IOringparam[3];	// �e��rms(�V�K�ǉ�)
//		//this->m  = Numeric::pi/6*pcd*pcd*pcd*FI.rhop[1];//�d�ʎ����v�Z
//		this->E  = IOringparam[1];
//		this->nu = IOringparam[2];
//		this->m  = Numeric::pi/6*pcd*pcd*pcd*IOringparam[0];	// �����O�����ǋ��̂ŉ���????
//
//	//case B4P_Ring::Outer:
//		//this->pcd  = (IOringparam[8] + IOringparam[9]) / 2 * 2;
//
//		//this->GV[0].x = -(IOringparam[6]+IOringparam[7])/2;
//		//this->GV[1].x = (IOringparam[6]+IOringparam[7])/2;
//		//this->GV[0].r = IOringparam[4];
//		//this->GV[1].r = IOringparam[5];
//		//this->GV[0].R = IOringparam[8];
//		//this->GV[1].R = IOringparam[9];
//		//this->GV[0].h = IOringparam[10];
//		//this->GV[1].h = IOringparam[11];
//		////??  = IOringparam[3];//�e��rms
//		////this->m  = Numeric::pi/6*pcd*pcd*pcd*FI.rhop[1];//�d�ʎ����v�Z
//		//this->E  = IOringparam[1];
//		//this->nu = IOringparam[2];
//
//		//this->m  = Numeric::pi/6*pcd*pcd*pcd*IOringparam[0];	// �����O�����ǋ��̂ŉ���D
//
//	//default:
//	//	break;
//	//}
//
//	double I_= this->m*pcd*pcd/10;
//	this->I << I_, I_, I_;
//	this->m_inv = 1.0 / this->m;
//	this->I_inv << 1 / I_, 1 / I_, 1 / I_;
//}


//// �������W�n���x���烊���O�̍a�������W�n���x�֕ϊ����郁�\�b�h�D
//double B4P_Ring::ine_to_Yvelocity
//	(					// out:	[m/s]:	�����O�̍a�������W�n�ɂ����鑬�x1�����D
//	const Vector3d& v,	// in : [m/s]:	�������W�n�̕ϊ����������x�D
//	double th			// in:	[rad]:	�����O�ʑ��p�D
//	)
//{
//	Vector3d u = this->get_cd(th);	// �ʑ��p�ł̐i�s�����x�N�g���D�i�������W�n�j
//	double   v_= u.dot(v);					// ���x�̐i�s���������D
//	return   v_;
//}
//
//// �������W�n���x���烊���O�̍a�������W�n���x�֕ϊ����郁�\�b�h�D
//Vector3d B4P_Ring::Y_to_inevelocity
//	(				// out:	[m/s]:	�������W�n�ɂ����鑬�x3�����D
//	double v,		// in : [m/s]:	�����O�̍a�������W�n�̕ϊ����������x�D
//	double th		// in:	[rad]:	�����O�ʑ��p�D
//	)
//{
//	Vector3d u  = this->get_cd(th);
//	Vector3d v_ = u * v;
//	return v_;
//}
//

//
//
//// ����ʑ��p�ł̎������x�N�g��(CircumferentialDirection)���擾���郁�\�b�h�D�iq=[1,0,0,0]�C �� = 45�� �� [0,-1/��2,-1/��2] �ƂȂ�j
//Vector3d B4P_Ring::get_cd
//	(				// out:	[m/s]:	�������W�n�ɂ����鑬�x3�����D
//	double th		// in:	[rad]:	�����ʒu�Ō����Ƃ��̃����O�ʑ��p�D
//	)
//{
//	Vector3d u_  = Vector3d(0.0, -cos(th), -sin(th));	// �����O���W�n�ł̐i�s�����x�N�g���D
//	Vector3d u__ = this->to_inevelocity(u_);			// �������W�n�ւ̕ϊ��D
//	return u__;
//}
//
//// ����ʑ��p�ł̌a�����x�N�g��(RadialDirection)���擾���郁�\�b�h�D�iq=[1,0,0,0]�C �� = 45�� �� [0,-1/��2,1/��2] �ƂȂ�j
//Vector3d B4P_Ring::get_rd
//	(				// out:	[m/s]:	�������W�n�ɂ����鑬�x3�����D
//	double th		// in:	[rad]:	�����ʒu�Ō����Ƃ��̃����O�ʑ��p�D
//	)
//{
//	Vector3d u_  = Vector3d(0.0, -sin(th), cos(th));	// �O�֍��W�n�ł̐i�s�����x�N�g���D
//	Vector3d u__ = this->to_inevelocity(u_);			// �������W�n�ւ̕ϊ��D
//	return u__;
//}


