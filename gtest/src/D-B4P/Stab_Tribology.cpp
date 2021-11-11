#include "pch.h"



// �g���N�V�����W���̐��l�Ƃ��ď��0.2���o�́D
double Tribology::Stab_Traction::calc
(						// out:	[-]:	�g���N�V�����W���D
	double eta0,		// in:	[Pas]:	���̔S�x�D
	double p0,			// in:	[Pa]:	�ő�ʈ��D
	double v0,			// in:	[m/s]:	���葬�x�D
	double u0			// in:	[m/s]:	�]���葬�x�D
) {
	return 0.2;
}

// �N�[�������C�W���̐��l�Ƃ��ď��0.3���o�́D�������C���葬�x��0�ɋ߂��ꍇ�C0���o�́D
double Tribology::Stab_Coulomb::calc(double mu0, double us, double ur, double s) {
	if (abs(us) < 1e-15) {
		return 0.0;
	}
	return 0.3;
}

// 
double Tribology::Stab_RollingResist::calc(double Rx, double Ry, double bl_D, double a, double b, double w, double E, double um, double fratio, double eta, double alpha, double beta, double k, double xio)
{
	return 0.1;
}