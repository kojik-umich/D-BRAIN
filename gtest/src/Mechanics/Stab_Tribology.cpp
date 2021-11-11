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
double Tribology::Stab_RollingResist::calc
(						// out:	[Nm]:	�g���N�D
	double Rx,			// in:	[m]:	�]������������ȗ��D
	double a,			// in:	[m]:	�ڐG�ȉ~���a�D
	double b,			// in:	[m]:	�ڐG�ȉ~�Z�a�D
	double w,			// in:	[N]:	�ڐG�׏d�D
	double E,			// in:	[Pa]:	����0,1�̓��������O���D
	double um,			// in:	[m/s]:	���ϑ��x�D
	double eta,			// in:	[Pas]:	���̔S�x�D
	double alpha,		// in:	[1/Pa]:	���̈��͔S�x�W���D
	double beta,		// in:	[1/K]:	���̉��x�S�x�W���D
	double k			// in:	[W/m2K]:���̔M�`�����D
) {
	return 0.1;
}