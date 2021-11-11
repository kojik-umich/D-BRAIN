#include "BS_Circuit.h"


BS_Circuit::BS_Circuit(void) {
	this->BL = NULL;
	return;
}


BS_Circuit::~BS_Circuit(void) {
	if (this->BL != NULL)
		delete[] this->BL;
	return;
}

void BS_Circuit::allocate(int n) {

	this->nBL = n;
	this->BL = new Ball[n];

	return;
}

void BS_Circuit::link(BS_Cylinder * CY, int is) {

	this->CY = CY;
	this->iSP = is;

	return;
}

void BS_Circuit::init(const BS_In::Circuit & circuit) {
	bool v_const[3], w_const[3];
	for (int i = 0; i < 3; i++) {
		v_const[i] = false;
		w_const[i] = false;
	}
	for (int i = 0; i < this->nBL; i++)

		this->BL[i].init(
			circuit.ball[i].r * 2,
			circuit.ball[i].young,
			circuit.ball[i].poisson,
			circuit.ball[i].density,
			circuit.ball[i].rms
			, v_const, w_const
		);
	this->th0 = circuit.th0;
	this->th1 = circuit.th1;

	return;
}

// �z������Ă���Cylinder����C���z�����ꂽ�ʍ��W�i�������W�n�j���擾���郁�\�b�h�D
void BS_Circuit::get_x0s(
	Vector3d*xs	// out:	[m]		
) {
	this->CY->linspace(this->iSP, this->th0, this->th1, this->nBL, xs);
	return;
}


// �ʂ̗����i�s�����itnb�̂���t�j���擾���郁�\�b�h�D���܂葬�x�͑����Ȃ��ł��D�J��Ԃ��v�Z�Ɏg���ɂ�Cylinder�̃��\�b�h�𒼐ڎg���Ă��������D
Vector3d BS_Circuit::get_t(int ib) {
	Vector3d eta = this->CY->to_etacoord(this->iSP, this->BL[ib].x);
	double th = eta[0];
	Matrix3d xyz2eta = this->CY->get_xyz2eta(this->iSP, th);
	return this->CY->to_inertialvector(Vector3d::UnitX(), xyz2eta);
}

double BS_Circuit::get_nd(void) {
	return this->CY->get_nd(this->iSP);
}

double BS_Circuit::get_r(void) {
	return this->CY->get_r(this->iSP);
}


