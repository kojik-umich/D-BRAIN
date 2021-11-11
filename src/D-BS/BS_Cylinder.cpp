/*******************************************************************************
!								"BS_Screw.cpp"
!													2020/01/08	[Core-T]	���
!	�˂����ނ�\���I�u�W�F�N�g�D
!	������p�����ăi�b�g�E�V���t�g���쐬����D
!	�a�̏𐔕����� new �ō쐬����D
!
!*******************************************************************************/
#include "BS_Cylinder.h"




// ���W���������W�n����i�b�g�̍a�������W�n�֕ϊ�
Vector3d BS_Cylinder::to_etacoord
(						// out:	[rad],[m]:	0�����F�i�b�g�ʑ��p�D1�����Feta���W�D2�����Fzeta���W�D
	int i,				// in : [-]:		�ϊ���������ԍ��D
	const Vector3d&x	// in : [m]:		�ϊ��������������W�n�D
) {
	Vector3d xyz = this->to_mycoord(x);
	Vector3d eta = this->SP[i].to_eta(xyz);
	return eta;
}

// ���W���i�b�g�̍a�������W�n���犵�����W�n�֕ϊ�
Vector3d BS_Cylinder::to_inertialcoord
(						// out:	[m]:		�������W�n�̍��W�D
	int i,				// in : [-]:		�ϊ���������ԍ��D
	const Vector3d&eta	// in :	[rad],[m]:	�ϊ��������������W�n�D
) {
	Vector3d xyz = this->SP[i].to_xyz(eta);
	Vector3d x = this->to_inecoord(xyz);
	return x;
}

// �x�N�g�����������W�n����i�b�g�̍a�������W�n�֕ϊ�
Vector3d BS_Cylinder::to_etavelocity
(						// out:	[rad],[m]:	0�����F�i�b�g�ʑ��p�D1�����Feta���W�D2�����Fzeta���W�D
	const Vector3d&v,	// in : [m/s]:		�ϊ����������x�i�������W�n�j�D
	const Matrix3d&xyz2eta
) {
	Vector3d v_ = this->to_myvelocity(v);
	Vector3d v__ = xyz2eta * v_;
	return v__;
}

// �x�N�g�����i�b�g�̍a�������W�n���犵�����W�n�֕ϊ�
Vector3d BS_Cylinder::to_inertialvelocity
(						// out:	[m]:		�������W�n�̍��W�D
	const Vector3d&v,	// in : [m/s]:		�ϊ����������x�i�������W�n�j�D
	const Matrix3d&xyz2eta
) {
	Vector3d v_ = xyz2eta.transpose() * v;
	Vector3d v__ = this->to_inevelocity(v_);
	return v__;
}

// �����ʑ��p����ϊ��s����Z�b�g���郁�\�b�h�D
Matrix3d BS_Cylinder::get_xyz2eta
(
	int i,				// in:	[-]:	�����ԍ�
	double theta		// in:	[rad]:	���݂̗����ʑ��p
) {
	Matrix3d xyz2eta = this->SP[i].get_xyz2eta(theta);
	return xyz2eta;
}


// �������W�n�̃x�N�g���𗆐����W�n�ɕϊ����郁�\�b�h�D
Vector3d BS_Cylinder::to_etavector(
	const Vector3d&a,	// in:	[any]	�ϊ��������x�N�g���i�������W�n�j
	const Matrix3d&xyz2eta
) {
	Vector3d b = this->to_myvector(a);
	Vector3d c = xyz2eta * b;

	return c;
}

// �������W�n�x�N�g�����������W�n�ɕϊ����郁�\�b�h�D
Vector3d BS_Cylinder::to_inertialvector(
	const Vector3d&a,	// in:	[any]	�ϊ��������x�N�g���i�������W�n�j
	const Matrix3d&xyz2eta
) {
	Matrix3d eta2xyz = xyz2eta.transpose();
	Vector3d b = eta2xyz * a;
	Vector3d c = this->to_inevector(b);
	return c;
}

// �ڐG�ȉ~���g�[���X�ɉ����ăX���C�X���郁�\�b�h�D
void BS_Cylinder::calc_slice
(
	int is,				// in:	[-]:	�����ԍ�;
	int ig,				// in:	[-]:	�a�ԍ��D0 or 1;
	double th,			// in:	[rad]:	�����ʑ��p�D
	const Vector3d&p,	// in:	[m]:	�ڐG�_�ʒu�̒��S�D�i�������W�n�j
	const Vector3d&xai,	// in:	[m]:	�a�f�ʂɂ����闆���i�s�����x�N�g���D�i�������W�n�j
	double a,			// in:	[m]:	�ڐG�ȉ~���a�D
	double b,			// in:	[m]:	�ڐG�ȉ~�Z�a�D
	const Vector2d*xy,	// in:	[-]:	�X���C�X���W�� [-] xy�͒P�ʉ~�̓����̒l�D0����:a�����D1����:b�����D
	int n,				// in:	[-]:	�ڐG�ȉ~�̃X���C�X���D
	const Vector2d&rho,	// in:	[1/m]:	�e�����ȗ��D0����:a�����D1����:b�����D
	Vector3d*ps			// out:	[m]:	�X���C�X���ꂽ�e�ڐG�_�ʒu�D�z��T�C�Y=n�D�i�������W�n�j
) {
	// �������W�n�ł́C�_p�����g�[���X�̐����f�ʉ~�̒��SO�D�i�g�[���X���o�i�i�Ɍ����Ă�ƁC�o�i�i�̐c�̃C���[�W�D�j
	Vector3d O = this->to_inertialcoord(is,
		Vector3d(th, this->SP[is].GV[ig].eta[0], this->SP[is].GV[ig].eta[1]));

	// �~�̌v�Z�̓R�X�g�������邽�߁C�ڐG�ȉ~���\�����������Ƃ���2���ߎ����s���D
	Vector3d op = p - O;
	double op_norm = op.norm();

	Vector3d cn = op / op_norm; // �a�f�ʂɐ����ȃx�N�g��
	Vector3d ct = cn.cross(xai);	// �ڐG�_���猩���aR�̎������D
	
	// �e�X���C�X�̒��S�_ pi �̍��W���Z�o
	for (int i = 0; i < n; i++) {
		double xi_ = a * xy[i].x();
		double yi_ = b * xy[i].y();
		double zi_ = - 0.5 * xi_ * xi_ * rho.x() - 0.5 * yi_ * yi_ * rho.y();
		ps[i] = p + cn * zi_ + ct * xi_ + xai * yi_;
	}
	return;
}

void BS_Cylinder::allocate(int ns) {

	this->nSP = ns;
	this->SP = new BS_Spiral[this->nSP];

	return;
};

// ���̓t�@�C�������Ƃɑ����������s���D
void BS_Cylinder::init(const BS_In::Cylinder & cylinder, const bool(&v_const)[3], const bool(&w_const)[3]) {

	this->set_const(v_const, w_const);

	this->rho = cylinder.density;
	this->nu = cylinder.poisson;
	this->E = cylinder.young;
	this->ri = cylinder.ri;
	this->ro = cylinder.ro;
	this->set_mI(cylinder.m, Vector3d(cylinder.Ix, cylinder.Iyz, cylinder.Iyz));

	for (int j = 0; j < this->nSP; j++)
		this->SP[j].init(cylinder.spiral[j]);

	return;
}

double BS_Cylinder::get_nd(int is) {
	return this->SP[is].get_nd();
}

double BS_Cylinder::get_r(int is) {
	return this->SP[is].get_r();
}

void BS_Cylinder::linspace
(
	int is,		// in:	[-]		��ԍ�
	double th0,	// in:	[rad]	�n�_�ʑ��p		
	double th1,	// in:	[rad]	�I�_�ʑ��p
	int nb,		// in:	[-]		�ʐ�
	Vector3d*xs	// out:	[m]		���z�����ꂽ�ʍ��W�i�������W�n�j
) {
	VectorXd th_ = VectorXd::LinSpaced(nb, th0, th1);

	for (int i = 0; i < nb; i++)
		xs[i] = this->to_inertialcoord(is, Vector3d(th_[i], 0.0, 0.0));

	return;
}

void BS_Cylinder::save(BS_Out::Cylinder & OUT) {

	Vector3d ax = this->get_ax();
	for (int k = 0; k < 3; k++) {
		OUT.x[k] = this->x[k];
		OUT.v[k] = this->v[k];
		OUT.w[k] = this->w[k];
		OUT.ax[k] = ax[k];
		OUT.F[k] = this->F[k];
		OUT.T[k] = this->T[k];
		OUT.Fs[k] = this->Fs[k];
		OUT.Ts[k] = this->Ts[k];
	}
	Vector4d q = this->q.coeffs();
	for (int k = 0; k < 4; k++)
		OUT.q[k] = q[k];

	return;
}


BS_Cylinder::BS_Cylinder() {
	this->SP = NULL;
	return;
}

BS_Cylinder::~BS_Cylinder() {
	if (this->SP != NULL)
		delete[] this->SP;
	return;
}

