/*******************************************************************************
!								"Spiral.cpp"
!													2019/12/25	[Core-T]	���
!	�����Ɋւ���v�Z�̃I�u�W�F�N�g�D
!	�������W�n�̗������W�n�̕ϊ�����Ȗ����D
!	����x���������Ă�����̂Ƃ���D�܂�x=0�ł̕��ʂł̎n�_�p�����Ƃ��C�E�˂��ŗ����͐i�ނ��̂Ƃ���D
!*******************************************************************************/
#include "Spiral.h"


void Spiral::init(
	double alp,		// in:	[rad]:	�����n�_�ʑ��p�D
	double l, 		// in:	[m]:	���[�h�����D
	double r		// in:	[m]:	�������a�D
) {
	this->alp = alp;
	this->l   = l;
	this->r   = r;

	this->l_2pi = l / (2 * Numeric::pi);
	this->nd    = sqrt(this->r*this->r + this->l_2pi*this->l_2pi);
	this->l_nd  = this->l_2pi / this->nd;
	this->r_nd  = this->r     / this->nd;

	return;
}

// ���ލ��W�n���痆�����W�n�ɕϊ����郁�\�b�h�D
Vector3d Spiral::to_eta
(						// out:	[rad],[m]	�������W�n
	const Vector3d&xyz	// in:	[m]:		���ލ��W�n�i����x�����������Ă��܂��j
) {
	double x      = xyz[0];
	double y      = xyz[1];
	double z      = xyz[2];
	double phi    = x / this->nd / this->l_nd;
	double beta   = phi + this->alp;
	double s_beta = sin(beta);
	double c_beta = cos(beta);
	double y0     = c_beta * y + s_beta * z;
	double z0     =-s_beta * y + c_beta * z;
	double th     = z0 / (y0 + this->l_nd * this->l_2pi / this->r_nd);
	double ze     = -th * this->l_2pi / this->r_nd;
	double et     = this->r_nd * this->nd - y0 - y0 * th*th / 2 + this->l_nd * th * ze;
	Vector3d eta  = Vector3d(th+phi, et, ze);
	return eta;
}

// �������W�n���畔�ލ��W�n�ɕϊ����郁�\�b�h�D
Vector3d Spiral::to_xyz
(							// out:	[m]:		�������W�n(x,y,z)(���̍��W�n)
	const Vector3d&eta		// in:	[rad],[m]	�������W�n(��,��,��)(���̍��W�n)
) {
	double th     = eta[0];
	double beta   = th + this->alp;
	double s_beta = sin(beta);
	double c_beta = cos(beta);
	Vector3d q    = Vector3d(this->l_2pi * th, this->r * c_beta, this->r * s_beta);
	Vector3d n    = Vector3d(0.0, -c_beta, -s_beta);
	Vector3d b    = Vector3d(this->r_nd, this->l_nd * s_beta, -this->l_nd * c_beta);
	Vector3d xyz  = q + eta[1] * n + eta[2] * b;
	return xyz;
}

// ���ލ��W�n���痆�����W�n�ɕϊ����郁�\�b�h�D1��J��Ԃ��v�Z���s�����ƂŁC���������̐��x�ŕϊ��ł���D
Vector3d Spiral::to_eta2
(							// out:	[rad],[m]:	�������W�n(��,��,��)(���̍��W�n)
	const Vector3d&xyz		// in:	[m]:		�������W�n(x,y,z)(���̍��W�n)
) {
	Vector3d eta0 = this->to_eta(xyz);
	Vector3d xyz0 = this->to_xyz(eta0);
	Vector3d xyz1 = 2 * xyz - xyz0;
	Vector3d eta1 = this->to_eta(xyz1);
	return eta1;
}

// x-y-z���W�n�̃x�N�g����xi-eta-zeta���W�n�̃x�N�g���ɕϊ����邽�߂̍s��擾�D
Matrix3d Spiral::get_xyz2eta
(							// out: [-]:	3*3�s��
	double theta			// in:	[rad]:	��ƂȂ闆�����W�n�̌��_�̈ʑ��p
) {
	double beta   = theta + this->alp;
	double s_beta = sin(beta);
	double c_beta = cos(beta);

	Matrix3d xyz2eta;
	xyz2eta <<
		this->l_nd, -this->r_nd * s_beta, this->r_nd * c_beta,
		0.0, -c_beta, -s_beta,
		this->r_nd, this->l_nd * s_beta, -this->l_nd * c_beta;
	return xyz2eta;
}

double Spiral::get_nd(void) {
	return this->nd;
}

double Spiral::get_r(void) {
	return this->r;
}
