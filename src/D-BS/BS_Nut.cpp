/*******************************************************************************
!								"BS_Nut.cpp"
!													2019/03/22	[Core-T]	���
!
!
!*******************************************************************************/
#include "BS_Nut.h"
using Eigen::VectorXd;


BS_Nut::BS_Nut() {
	this->CY = NULL;
	return;
}

BS_Nut::~BS_Nut() {
	if (this->CY != NULL)
		delete[] this->CY;
	return;
}

// �ȉ��C�V���O���i�b�g�D
void BS_SingleNut::Nothing(void) {
	return;
}

void BS_Nut::allocate(const std::vector<BS_In::Cylinder>&cylinders) {

	this->nCY = cylinders.size();
	this->CY = new BS_Cylinder[this->nCY];

	for (int i = 0; i < this->nCY; i++)
		this->CY[i].allocate(cylinders[i].spiral.size());
	return;
}

void BS_SingleNut::init(const std::vector<BS_In::Cylinder> & nut, double w0) {
	bool v_const[3], w_const[3];
	for (int i = 0; i < 3; i++) {
		v_const[i] = true;
		w_const[i] = true;
	}

	this->set_const(v_const, w_const);


	Vector3d x, v, w; Quaterniond q;
	x = Vector3d::Zero();
	v = Vector3d::Zero();
	w = Vector3d::Zero();
	q = Quaterniond::Identity();

	this->set_y_(x, v, q, w);

	this->CY[0].init(nut[0], v_const, w_const);
	return;
}

void BS_SingleNut::set_y_(const Vector3d & x, const Vector3d & v, const Quaterniond & q, const Vector3d & w) {
	this->set_y(x, v, q, w);
	this->CY[0].set_y(x, v, q, w);
	return;
}

// �ȉ��C�_�u���i�b�g�D
void BS_DoubleNut::Nothing(void) {
	return;
}



void BS_DoubleNut::init(const std::vector<BS_In::Cylinder> & nut, double w0) {
	bool v_const[3], w_const[3];
	for (int i = 0; i < 3; i++) {
		v_const[i] = true;
		w_const[i] = true;
	}
	this->set_const(v_const, w_const);

	for (int i = 0; i < 2; i++)
		this->dx[i] = Vector3d(nut[i].x0, 0.0, 0.0);

	Vector3d x, v, w; Quaterniond q;
	x = Vector3d::Zero();
	v = Vector3d::Zero();
	w = Vector3d(w0, 0.0, 0.0);
	q = Quaterniond::Identity();

	this->set_y(x, v, q, w);
	this->set_dx();

	for (int i = 0; i < 2; i++)
		this->CY[i].init(nut[i], v_const, w_const);

	return;
}

// �i�b�g���S���獶�E�̃i�b�g���ꂼ��̕ψʂ�^���郁�\�b�h�D
void BS_DoubleNut::set_dx(void) {

	for (int i = 0; i < 2; i++) {
		Vector3d dx = this->to_inecoord(this->dx[i]);
		this->CY[i].set_param(this->x + dx, this->v, this->q, this->w);
	}
	return;
}

// �l�̐ݒ胁�\�b�h�D�������ʂœn���ꂽ�l��L�����ɂ��ăp�����^�Ƃ��đ������D
void BS_DoubleNut::set_y_(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w) {

	this->set_y(x, v, q, w);
	this->set_dx();

	return;
}






























//#include "BS_Screw.h"
//
//
//void BS_Screw::init(int ns) {
//
//	this->SP = new BS_Spiral[ns];
//}
//
//// �i�b�g�̍a�������W�n�ɕϊ����郁�\�b�h�D
//Vector3d BS_Screw::to_etacoord
//(						// out:	[rad],[m]:	0�����F�i�b�g�ʑ��p�D1�����Feta���W�D2�����Fzeta���W�D
//	const Vector3d&x,	// in : [m]:		�ϊ��������������W�n�D
//	int i				// in : [-]:		�ϊ���������ԍ��D
//) {
//	Vector3d xyz = this->to_mycoord(x);
//	Vector3d eta = this->SP[i].to_eta2(xyz);
//	return eta;
//}
//
//// �i�b�g�̍a�������W�n�ɕϊ����郁�\�b�h�D
//Vector3d BS_Screw::to_inertialcoord
//(						// out:	[m]:		�������W�n�̍��W�D
//	const Vector3d&eta,	// in :	[rad],[m]:	�ϊ��������������W�n�D
//	int i				// in : [-]:		�ϊ���������ԍ��D
//) {
//	Vector3d xyz = this->SP[i].to_xyz(eta);
//	Vector3d x = this->to_inecoord(xyz);
//	return x;
//}
// 
//void BS_Screw::set_TransMat
//(
//	double theta,	// in : [rad]:		�i�b�g�ʑ��p�D
//	int i			// in : [-]:		�ϊ���������ԍ��D
//)
//{
//	this->SP[i].set_TransMat(theta);
//}
//// �������W�n���x���痆�����W�n���x�֕ϊ����郁�\�b�h�D
//// ���g���O�ɁC�K��set_TransMat�ōX�V���Ă��������I
//Vector3d BS_Screw::to_etavelocity
//(					// out:	[m/s]:	�������W�n�ɂ����鑬�x3�����D
//	Vector3d v,		// in : [m/s]:	�������W�n�̕ϊ����������x�D
//	int i			// in : [-]:	�ϊ���������ԍ��D
//)
//{
//	Vector3d v_ = this->to_myvelocity(v);			// �i�b�g���W�n�ւ̕ϊ��D
//	Vector3d v__ = this->SP[i].to_etavector(v_);	// �������W�n�ւ̕ϊ��D
//	return v__;
//}
//
//// �������W�n���x���痆�����W�n���x�֕ϊ����郁�\�b�h�D
//// ���g���O�ɁC�K��set_TransMat�ōX�V���Ă��������I
//Vector3d BS_Screw::to_inertialvelocity
//(					// out:	[m/s]:	�������W�n�ɂ����鑬�x3�����D
//	Vector3d v,		// in : [m/s]:	�������W�n�̕ϊ����������x�D
//	int i			// in : [-]:	�ϊ���������ԍ��D
//)
//{
//	Vector3d v_ = this->SP[i].to_xyzvector(v);		// �i�b�g���W�n�ւ̕ϊ��D
//	Vector3d v__ = this->to_inevelocity(v_);		// �������W�n�ւ̕ϊ��D
//	return v__;
//}

