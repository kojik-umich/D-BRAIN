#include "Ball.h"

// �_�~�[�֐��D
void Ball::Nothing(void) {
}


// Ball([��]0�ʔԍ��C1���x�C2�����O���C3por�C4�e��rms�C5�ʌa)�̃p�����[�^����
void Ball::init(
	double D,		// �ʌa
	double E,		// �����O��
	double por,		// �|�A�\����
	double den,		// ���x
	double rms,		// �e��rms
	const bool(&v_const)[3],	// �S�������i���x�j
	const bool(&w_const)[3]		// �S�������i�p���x�j
) {
	this->D				= D;	
	this->E				= E;	
	this->nu			= por;	
	this->rho			= den;	
	this->sigmap		= rms;	
	this->m				= Numeric::pi / 6 * D * D * D * rho;
	this->m_inv			= 1.0 / this->m;
	this->I				= this->m*D*D/10 * Vector3d::Ones();
	this->I_inv			= 1.0 / this->I[0] * Vector3d::Ones();
	this->r				= 0.5 * D;
	this->r_inv			= 2.0 / D;
	this->set_const(v_const, w_const);
	return;
}

// �ʂɐڋ߂��镨�̂̍��W�E���x�����ɁC�ڐG�̗L���C�ڋߗʁC�ڋߑ��x�C�ڋߕ������v�Z
bool Ball::calc_Contact(	// out:	[-]:	�ڐG�̗L���i�ڐG����Ftrue�C�ڐG�Ȃ��Ffalse�j
	const Vector3d&x,		// in:	[m]:	�ڋ߂��镨�̂̍��W�D�i�������W�n�j
	const Vector3d&v,		// in:	[m/s]:	�ڋ߂��镨�̂̕\�ʑ��x�D�i�������W�n�j
	double &dx,				// out:	[m]:	�ڋߗʁi�X�J���[�C�����Ă���ΐ��j
	double &dv,				// out:	[m/s]:	�ڋߑ��x�i�X�J���[�C�ڋߗʂ��傫���Ȃ���������j
	Vector3d&edir			// out:	[-]:	�ڋߕ����P�ʃx�N�g���D
) {

	Vector3d Zero = Vector3d::Zero();

	// �{�[�����S����N���҂ւ̋����ƁC���̕����̒P�ʃx�N�g���D
	double Distance = (x - this->x).norm();
	if (Distance == 0)
		return false;

	// �{�[�����S����N���҂Ɍ������x�N�g���D�{�[���ɂ�����׏d�Ƌt�����D
	edir = (x - this->x) / Distance;

	dx = this->r - Distance;

	// �ڐG���Ă��Ȃ��ꍇ�׏d0�ŕԂ��D
	if (dx < 0)
		return false;

	// �ڋ߂��镨�̂̑��x dv ���Z�o����D�i�X�J���[�C�ڋ߂�����������j
	dv = (v - this->v).dot(-edir);

	return true;
}



// ���͂��ꂽ�x�N�g���́C�@�������������������郁�\�b�h�D
Vector3d Ball::remove_Normal(	// out:	[]:		�@�����������������x�N�g���D�i�������W�n�j
	const Vector3d&p,			// in:	[m]:	�N���҂̍��W�D�i�������W�n�j
	const Vector3d&a			// in:	[]:		���H�������x�N�g���D�i�������W�n�j
) {
	Vector3d e = (p - this->x).normalized();
	double ea = e.dot(a);
	return a - ea * e;
}

void BallBallPair::link(Ball * BL0, Ball * BL1) {
	this->BL[0] = BL0;
	this->BL[1] = BL1;
	return;
}
