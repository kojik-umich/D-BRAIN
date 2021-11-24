/*******************************************************************************
!								"Tribology.cpp"
!													2018/11/27	[Core-T]	���
!
!
!
!
!
!*******************************************************************************/

#include "Tribology.h"

// �ڐG�׏d�̃X���C�X������ϕ��ɂ��v�Z�D
void Tribology::SliceForceRatio
(
	int msmax,			// in:	[-]:	�ڐG�ȉ~���a�����ɃX���C�X���鐔�D
	double*ForceRatio	// out:	[N]:	�ڐG�׏d�������z�D�z��T�C�Y��msmax�D
	) {
	for (int i = 0; i < msmax; i++) {
		double x1 = -1 + 2 * i / (double)msmax;					// i�Ԗڂ̗v�f��x���W��Z���̕��Ŋ��������l
		double x2 = -1 + 2 * (i + 1) / (double)msmax;
		ForceRatio[i] = 0.25 * (3 * (x2 - x1) - (Numeric::Cube(x2) - Numeric::Cube(x1)));
	}
	return;
}

// �ڐG�׏d�̃X���C�X������ϕ��ɂ��v�Z�D
void Tribology::SliceForceRatio2d
(
	int nr,		// in:	[-]:	�ڐG�ȉ~r�����ɃX���C�X���鐔�D
	int nt,		// in:	[-]:	�ڐG�ȉ~�ƕ����ɃX���C�X���鐔�D
	double*r	// out:	[N]:	�ڐG�׏d�������z�D�z��T�C�Y:nr�D
) {
	double dt = 2 * Numeric::pi / nt;

	for (int i = 0; i < nr; i++) {
		double x0 = double(i)     / nr;
		double r0 = -pow((1 - x0 * x0), 1.5) / 2 / Numeric::pi;

		double x1 = double(i + 1) / nr;
		double r1 = -pow((1 - x1 * x1), 1.5) / 2 / Numeric::pi;
		
		r[i] = (r1 - r0) * dt;
	}
	return;
}

// ���������O�����Z�o����D
double Tribology::ReducedYoung
(						// out:	[Pa]:	���������O���D
	double E0,			// in:	[Pa]:	����0�̃����O���D
	double m0,			// in:	[-]:	����0�̃|�A�\����D
	double E1,			// in:	[Pa]:	����1�̃����O���D 
	double m1			// in:	[-]:	����0�̃|�A�\����D
) {
	double E = 2.0 / ((1.0 - m0*m0)/E0 + (1.0 - m1*m1)/E1);
	return E;
}

// ���Z���ʂ��Z�o����D
double Tribology::ReducedMass
(						// out:	[kg]:	���Z���ʁD
	double m0,			// in:	[kg]:	����0�̎��ʁD
	double m1			// in:	[kg]:	����1�̎��ʁD
) {
	double M = 1.0 / (1.0 / m0 + 1.0 / m1);
	return M;
}

// �����e�����Z�o����D
double Tribology::CompositeRoughness
(						// out:	[m]:	�����e���D
	double s0,			// in:	[m]:	����0�̑e���D
	double s1			// in:	[m]:	����1�̑e���D
) {
	double s = sqrt(s0 * s0 + s1 * s1);
	return s;
}

// �����_�l�ɂ���Ė����׏d�������Z�o����D
double Tribology::ForceRatio
(						// out:	[-]:	�����׏d�x�������D
	double lambda		// in:	[-]:	�����_�l�D
) {
	if (lambda < 0)
		return 0.0;

	return 1.0 - exp(-1.8 * pow(lambda, 1.2));
}
// ��-�a�ڐG���̐ڐG�ȉ~���a���炱��L�������i�̔����j���v�Z�i������457���j
double Tribology::RollerLength(
	double a,	// �ڐG�������a [m]
	double Rx,	// �]��������������a [m]
	double Ry,	// �]��������ƒ�����������̓������a [m]
	double xio	// ���j�X�J�X���� [m]
) {

	// ���j�X�J�X����������`�̂Ƃ���a��1.2�{�����뒷���Ƃ���D
	if (xio < 0) {
		double wi = 1.2 * a;
		return wi;
	}
	double Xio = xio / Rx;
	double A = a / Rx;
	double Hio = 0.5 * Xio * Xio;
	double ayw = Ry / Rx;

	double wi = Ry * sqrt(2 * ayw * Hio + A * A);
	// ����L��������a��1.2�{�𒴂��Ȃ��悤�ɂ���i���R�͌����Q�Ɓj
	if (wi > 1.2 * a) {
		wi = 1.2 * a;
	}

	return wi;
}

// �^���W�F���g�J�[�u�̖��C�W���D
double Tribology::Tangent::calc(double mu0, double us, double ur, double s) {

	return mu0 * 2.0 / Numeric::pi * atan(s * us / ur);
}

// �N�[�������C���Ȃ���Ԃ��`�D
double Tribology::CoulombNothing::calc(double mu0, double us, double ur, double s) {
	return 0.0;
}


double Tribology::SimpleCoulomb::calc(double mu0, double us, double ur, double s) {
	return mu0;
}

// Brew Hamrock �̋ߎ��v�Z�D���ȏ��̒l�Ƃ̈�v���m�F�D
void Tribology::BrewHamrock::calc
(
	double Rx,			// in:	[m]:	�����������ȗ��D
	double Ry,			// in:	[m]:	�]������������ȗ��D
	double dx,			// in:	[m]:	�{�[���H�����ݗʁD
	double E,			// in:	[Pa]:	���������O���D
	double&k,			// out:	[N/m^1.5]:	����`�����D
	double&a,			// out:	[m]:	�ڐG�ȉ~���a�D
	double&b			// out:	[m]:	�ڐG�ȉ~�Z�a�D
) {
	double sum_rho = 1. / Rx  + 1. / Ry;
	double Kk = 1.5277 + 0.6023 * log(Ry/Rx);										// ���튮�S�ϕ�
	double Ek = 1.0003 + 0.5968 * Rx / Ry;											// ���튮�S�ϕ�
	double mu = 0.87959 * pow(Ek, Numeric::ot) * pow((Ry/Rx), 0.4240);				// �W�� ��

	double k_el = 1 / (1.0339 * pow(Ry/ Rx, 0.636));								// �ȉ~�� k_el 
	double nu = pow(2 * Ek * k_el / Numeric::pi, Numeric::ot);						// �W�� ��

	double K2pimu = 2 * Kk / Numeric::pi / mu;										// �W�� K/2�΃�

	k = pow(pow(1.125 / E / E * sum_rho, Numeric::ot) * K2pimu, -1.5);				// �e���W���̎Z�o
	a = pow(6 * Ek * k / E / sum_rho / Numeric::pi / k_el / k_el , Numeric::ot) * sqrt(dx);	// �����a a �̎Z�o
	b = k_el * a;																	// �Z���a b �̎Z�o
	return;
}

// Hydrodynamic force and moment in pure rolling lubricated contacts. Part 2: Point contacts �Ɋ�Â��C�g���N�v�Z���s���D
double Tribology::Houpert::calc
(						// out:	[Nm]:	�g���N�D
	double Rx,			// in:	[m]:	�]������������ȗ��D
	double Ry,			// in:	[m]:	�a���������ȗ��D
	double bl_D,		// in:	[m]:	�ʌa
	double a,			// in:	[m]:	�ڐG�ȉ~���a�D
	double b,			// in:	[m]:	�ڐG�ȉ~�Z�a�D
	double w,			// in:	[N]:	�ڐG�׏d�D
	double E,			// in:	[Pa]:	����0,1�̓��������O���D
	double um,			// in:	[m/s]:	���ϑ��x�D
	double fratio,		// in:	[-]:	�����׏d�����D
	double eta,			// in:	[Pas]:	���̔S�x�D
	double alpha,		// in:	[1/Pa]:	���̈��͔S�x�W���D
	double beta,		// in:	[1/K]:	���̉��x�S�x�W���D
	double k,			// in:	[W/m2K]:���̔M�`�����D
	double xio			// in:	[m]:	���j�X�J�X����
) {
	double p     = 1.5 * w / (Numeric::pi * a * b);
	double kap   = b / a;
	double beta_ = b * E / (Numeric::pi * p * Rx);
	double W     = 2 * Numeric::pi * p * b * b / (3 * kap * E * Rx * Rx);
	double U     = 2 * eta * um / E / Rx;
	double M     = pow(beta_, 0.25) * kap * W * pow(U, -0.75);
	double T     = 2 * a / Rx * ((0.77 * pow(beta_, -1/3) * pow(kap, 0.12) * pow(U, -1./12.) -1.8 * pow(beta_, -0.25)) / (1 + M / 6.6) + 1.8 * pow(beta_, -0.25)) * pow(U, 0.75);

	double t     = T * E * Rx * Rx * Rx;
	return t;
}

// EHL rolling resistance formula in Aihara's paper [ASME J. Trib., 109 (1987), 471] �̎�
double Tribology::AiharaR::calc
(						// out:	[Nm]:	�g���N�D
	double Rx,			// in:	[m]:	�]������������ȗ��D
	double Ry,			// in:	[m]:	�a���������ȗ��D
	double bl_D,		// in:	[m]:	�ʌa
	double a,			// in:	[m]:	�ڐG�ȉ~���a�D
	double b,			// in:	[m]:	�ڐG�ȉ~�Z�a�D
	double w,			// in:	[N]:	�ڐG�׏d�D
	double E,			// in:	[Pa]:	����0,1�̓��������O���D
	double um,			// in:	[m/s]:	���ϑ��x�D
	double fratio,		// in:	[-]:	�����׏d�����D
	double eta,			// in:	[Pas]:	���̔S�x�D
	double alpha,		// in:	[1/Pa]:	���̈��͔S�x�W���D
	double beta,		// in:	[1/K]:	���̉��x�S�x�W���D
	double k,			// in:	[W/m2K]:���̔M�`�����D
	double xio			// in:	[m]:	���j�X�J�X����
) {
	double U = eta * um / E / Rx;
	double G = alpha * E;
	double W = w / E / Rx / Rx;
	//double W = 0.75 * w / E / Rx / a;										// ���ڐG�̉׏d�p�����[�^�̋��ߕ��ɍ��킹��
	//double lc = a;// Tribology::RollerLength(a, Rx, Ry, xio);				// ����̗L������
	double L = eta * beta * um * um / k;
	double M = 176 / (1 + 0.29 * pow(L, 0.78)) / alpha * pow(G * U, 0.658) * pow(W, 0.31) * Rx * Rx * a;
	//double M = 176 / (1 + 0.29 * pow(L, 0.78)) / alpha * pow(G * U, 0.658) * pow(W, 0.31) * Rx * Rx * lc;
	return M;
}

// NTN TECHNICAL REVIEW No.82�i2014�j���C�����̎��D
double Tribology::Fujiwara::calc
(						// out:	[Nm]:	�g���N�D
	double Rx,			// in:	[m]:	�]������������ȗ��D
	double Ry,			// in:	[m]:	�a���������ȗ��D
	double bl_D,		// in:	[m]:	�ʌa
	double a,			// in:	[m]:	�ڐG�ȉ~���a�D
	double b,			// in:	[m]:	�ڐG�ȉ~�Z�a�D
	double w,			// in:	[N]:	�ڐG�׏d�D
	double E,			// in:	[Pa]:	����0,1�̓��������O���D
	double um,			// in:	[m/s]:	���ϑ��x�D
	double fratio,		// in:	[-]:	�����׏d�����D
	double eta,			// in:	[Pas]:	���̔S�x�D
	double alpha,		// in:	[1/Pa]:	���̈��͔S�x�W���D
	double beta,		// in:	[1/K]:	���̉��x�S�x�W���D
	double k,			// in:	[W/m2K]:���̔M�`�����D
	double xio			// in:	[m]:	���j�X�J�X����
) {
	double U   = eta * um / E / Rx;
	double G   = alpha * E;
	double W   = w / E / Rx / Rx;
	double kap = b / a;
	double FHD = 44.6 * pow(U, 0.694) * pow(G, 0.961) * pow(W, 0.401) * (1 - 0.962 * exp(-0.818 * kap));
	double Fr  = FHD * Rx * Rx / alpha;

	double M  = Fr * Rx;
	return M;
}

// Baltac�Ɏ�������Ă���"�����̎�"�D�_���f�ڂ̑����̎��ƂƓ��l��Goksem�̎������ǂ��č쐬�������̂����C�ʎ���p�Ƀ`���[�j���O����Ă���i�����񍐑�457���ɋL�ځj
double Tribology::GoksemAihara::calc
(						// out:	[Nm]:	�g���N�D
	double Rx,			// in:	[m]:	�]������������ȗ��D
	double Ry,			// in:	[m]:	�a���������ȗ��D
	double bl_D,		// in:	[m]:	�ʌa
	double a,			// in:	[m]:	�ڐG�ȉ~���a�D
	double b,			// in:	[m]:	�ڐG�ȉ~�Z�a�D
	double w,			// in:	[N]:	�ڐG�׏d�D
	double E,			// in:	[Pa]:	����0,1�̓��������O���D
	double um,			// in:	[m/s]:	���ϑ��x�D
	double fratio,		// in:	[-]:	�����׏d�����D
	double eta,			// in:	[Pas]:	���̔S�x�D
	double alpha,		// in:	[1/Pa]:	���̈��͔S�x�W���D
	double beta,		// in:	[1/K]:	���̉��x�S�x�W���D
	double k,			// in:	[W/m2K]:���̔M�`�����D
	double xio			// in:	[m]:	���j�X�J�X����
) {

	double lc =  Tribology::RollerLength(a, Rx, Ry, xio);		// ����̗L�������i�̔���?�j
					
	double w_mj = 0.75 * w / a;			// ���׏d[N/m]
	double Wehl = w_mj / (E * Rx);		// �׏d�p�����[�^[-]
	double G = alpha * E;
	double U = eta * um / E / Rx;

	double L_GH = eta * beta * um * um / k;			// �M���׌W��
	double D_GH = 11.81220746 * G * U / pow(Wehl, 1.5);			// �����׌W��
	double Kt = exp((0.658228 + 0.0106383 * pow(L_GH, 0.3174603)) * log(D_GH) - (0.162625 + log(1 + 0.411483 * pow(L_GH, 0.707071))));
	
	double el = 2 * lc * 1000;					// ����̗L������[mm]
	double Da = bl_D * 1000;				// �ʌa[mm]
	double wle = w / 9.8065 / el;
	double f_Q = pow((2 * fratio * wle / Da), 0.3);	// �׏d�␳���i���������ɂȂ�͂������C�Ȃ����������ł͂Ȃ��j

	
	double Cst = this->Starvation(xio, D_GH, L_GH, b, Kt);	// �X�^�x�[�V�����W��

	double tr = Cst * (1 / alpha) * fratio * Wehl * Kt * Rx;
	double Mr = tr * 2 * Rx * 2 * lc * f_Q;

	return Mr;
}

// �q�X�e���V�X�������v�Z�D�v�Z���� "�p�c�a�Y, �ʎ���̖��C���[�����g�Ɋւ��錤���i�X���X�g����̏ꍇ�j�@�B�w�27��178��"�Ɍf��
double Tribology::Kakuta::calc
(					// out: [Nm]:	�q�X�e���V�X����
	double fh,		// in:	[-]:	�q�X�e���V�X�W��
	double b,		// in:	[m]:	�ڐG�ȉ~�Z���a
	double Q		// in:	[N]:	�]���̉׏d
) {
	double Mrh = 0.375 * fh * b* Q;
	return Mrh;
}

// �q�X�e���V�X�������Ȃ��Ƃ�
double Tribology::HysteresisNothing::calc
(					// out: [Nm]:	�q�X�e���V�X����
	double fh,		// in:	[-]:	�q�X�e���V�X�W��
	double b,		// in:	[m]:	�ڐG�ȉ~�Z���a
	double Q		// in:	[N]:	�]���̉׏d
) {
	return 0;
}

// Goksem�̎��̃X�^�x�[�V�����W��
double Tribology::GoksemAihara::Starvation(
	double xio,		//	in:	[m]:	���j�X�J�X����
	double D_GH,	//	in: [-]:	�����׌W��
	double L_GH,	//	in: [-]:	�M���׌W��
	double b,		//	in:	[m]:	�ڐG�ȉ~�Z���a
	double Kt		//	in:	[-]:	??
){
	// �\������ or �X�^�x�[�V�����W�����ړ��͂̂Ƃ��͕␳����
	if (xio <= 0) {
		return 1;
	}
	double Ca = 950 * pow(D_GH, 0.05) / (1000 + (0.2 * D_GH) * pow(L_GH, (0.66 - 0.003 * D_GH))); // �X�^�x�[�V�����̉e��[-]
	double Cb = pow(5 + D_GH, 0.068) * (pow(1 + L_GH, 0.06 / pow(D_GH, 0.3) - 0.035) - 0.05 * pow(L_GH, 1e-4));
	double mst = 1 + xio / b;
	double Kai = mst * sqrt(mst * mst - 1) - log(mst + sqrt(mst * mst - 1));
	double fai = pow(3 * Kai / (42 * Kt), 0.3);
	double Cst = Ca * pow(fai, Cb) / (1 + Ca * pow(fai, Cb));

	return Cst;

}


// �]���薀�C���Ȃ���Ԃ��`�D
double Tribology::RollingResistanceNothing::calc(double Rx, double Ry, double bl_D, double a, double b, double w, double E, double um, double fratio, double eta, double alpha, double beta, double k, double xio) {
	return 0.0;
}



// HamrockDowson�̒������������̎��D���M�ɂ��S�x�ω����l���������ꍇ�͉���ErtelGrubin���������Ă��������D�i���`�v�Z�Ȃ̂ł��̂܂܊|���Z�ł��܂��D�j�i�o�T�FHamrock, Bernard J., and D. Dowson. "Isothermal elastohydrodynamic lubrication of point contacts. IV-Starvation results." (1976).�j
double Tribology::HamrockDowsonHc::calc
(						// out:	[m]:	���������D
	double w,			// in:	[N]:	�ڐG�׏d�D
	double E,			// in:	[Pa]:	���������O���D
	double Rx,			// in:	[m]:	�]��������ȗ����a�D
	double Ry,			// in:	[m]:	�������ȗ����a�D
	double alpha,		// in:	[1/Pa]:	���͔S�x�W���D
	double eta,			// in:	[Pas]:	���̔S�x�D
	double u,			// in:	[m/s]:	�]���葬�x�D
	double lm,			// in:	[m]:	���j�X�J�X�����D�i<-1:�����v�Z�C1~0:���̂܂܊|���Z�C>0:���j�X�J�X�����g�p�D�j
	double bh			// in:	[m]:	�ڐG�ȉ~�Z�a�D
) {
	double W = w / E / Rx / Rx;
	double G = alpha * E;
	double U = eta * u / E / Rx;
	double k_ = 1.0339 * pow(Ry / Rx, 0.636);

	double H = 2.69 * pow(U, 0.67) * pow(G, 0.53) * pow(W, -0.067) * (1.0 - 0.61 * exp(-0.73 * k_));
	double h =  H * Rx;

	h *=  this->Starvation(H, Rx, lm, bh);
	return h;
}

// HamrockDowson�̍ŏ����������̎��D���M�ɂ��S�x�ω����l���������ꍇ�͉���ErtelGrubin���������Ă��������D�i���`�v�Z�Ȃ̂ł��̂܂܊|���Z�ł��܂��D�j�i�o�T�FHamrock, Bernard J., and D. Dowson. "Isothermal elastohydrodynamic lubrication of point contacts. IV-Starvation results." (1976).�j
double Tribology::HamrockDowsonHmin::calc
(						// out:	[m]:	���������D
	double w,			// in:	[N]:	�ڐG�׏d�D
	double E,			// in:	[Pa]:	���������O���D
	double Rx,			// in:	[m]:	�]��������ȗ����a�D
	double Ry,			// in:	[m]:	�������ȗ����a�D
	double alpha,		// in:	[1/Pa]:	���͔S�x�W���D
	double eta,			// in:	[Pas]:	���̔S�x�D
	double u,			// in:	[m/s]:	�]���葬�x�D
	double lm,			// in:	[m]:	���j�X�J�X�����D�i<-1:�����v�Z�C1~0:���̂܂܊|���Z�C>0:���j�X�J�X�����g�p�D�j
	double bh			// in:	[m]:	�ڐG�ȉ~�Z�a�D
) {
	double W = w / E / Rx / Rx;
	double G = alpha * E;
	double U = eta * u / E / Rx;
	double k_ = 1.0339 * pow(Ry / Rx, 0.636);

	double H = 3.63 * pow(U, 0.68) * pow(G, 0.49) * pow(W, -0.073) * (1.0 - exp(-0.68 * k_));
	double h = H * Rx;

	h *= this->Starvation(H, Rx, lm, bh);
	return h;
}



// HamrockDowson�������������̃X�^�x�[�V�����W���D
double Tribology::HamrockDowsonHc::Starvation
(						// out:	[-]:	�X�^�x�[�V�����W���D
	double H,			// in:	[-]:	�����p�����^�D
	double Rx,			// in:	[m]:	�]��������ȗ����a�D
	double lm,			// in:	[m]:	���j�X�J�X�����D�i<-1:�����v�Z�C1~0:���̂܂܊|���Z�C>0:���j�X�J�X�����g�p�D�j
	double bh			// in:	[m]:	�ڐG�ȉ~�Z�a�D
) {
	if (lm <=  -1.0)
		return 1.0;

	if (lm <=  0.0)
		return fabs(lm);

	double f  = lm / bh;
	double fc = 3.06 * pow(H * Rx * Rx / bh / bh, 0.58);
	if (f > fc)
		return 1.0;

	return pow(f / fc, 0.29);
}

// HamrockDowson�ŏ����������̃X�^�x�[�V�����W���D
double Tribology::HamrockDowsonHmin::Starvation
(						// out:	[-]:	�X�^�x�[�V�����W���D
	double H,			// in:	[-]:	�����p�����^�D
	double Rx,			// in:	[m]:	�]��������ȗ����a�D
	double lm,			// in:	[m]:	���j�X�J�X�����D�i<-1:�����v�Z�C1~0:���̂܂܊|���Z�C>0:���j�X�J�X�����g�p�D�j
	double bh			// in:	[m]:	�ڐG�ȉ~�Z�a�D
) {
	if (lm <= -1.0)
		return 1.0;

	if (lm <= 0.0)
		return fabs(lm);

	double f = lm / bh;
	double fc = 3.34 * pow(H * Rx * Rx / bh / bh, 0.58);
	if (f > fc)
		return 1.0;

	return pow(f / fc, 0.25);
}

// �����������Ȃ���Ԃ��`�D
double Tribology::FilmThicknessNothing::calc(double w, double E, double Rx, double Ry, double alpha, double eta, double u, double lm, double bh) {
	return 0.0;
}

// ���M�ɂ��S�x�ω����l�������␳�W���D
double Tribology::ErtelGrubin
(						// out:	[-]:	�␳�W���i0~1�̎����j�D���������ɂ����Z���s���D
	double eta,			// in:	[Pas]:	���̔S�x�D
	double beta,		// in:	[1/K]:	���x�S�x�W���D
	double k,			// in:	[W/m2K]:���̔M�`�����D
	double u			// in:	[m/s]:	�]���葬�x�D
) {
	double L = eta * beta * u * u / k;
	double phi =  3.94 / (3.94 + pow(L, 0.62));
	return phi;
}

// �����̎�����g���N�V�����W�����v�Z����D
double Tribology::AiharaT::calc
(						// out:	[-]:	�g���N�V�����W���D
	double eta0,		// in:	[Pas]:	���̔S�x�D
	double p0,			// in:	[Pa]:	�ő�ʈ��D
	double v0,			// in:	[m/s]:	���葬�x�D
	double u0			// in:	[m/s]:	�]���葬�x�D
) {
	if (eta0 <= 0.0)
		return 0.0;

	if (p0 <= 0.0)
		return 0.0;

	if (v0 <= 0.0)
		return 0.0;

	if (u0 <= 0.0)
		return 0.0;

	double s    = v0 / u0;
	double eta  = Unit::Pas2cP(eta0);
	double u    = Unit::ms2mms(u0);
	double pmax = Unit::Pa2kgfmm2(p0);
	double smax = 1.2448e7 * pow(eta, -1.5) * pow(u, -0.75) * pow(pmax, -1.5);

	double mumax;
	if (pmax > 94.404)
		mumax = 4.078e-2 * pow(eta, 0.1) * pow(u, -0.12) * pow(pmax, 0.18);
	else
		mumax = 1.008e-4 * pow(eta, 0.1) * pow(u, -0.12) * pow(pmax, 1.5);

	double mu = mumax * (1.0 - exp(-2.6 * s / smax));

	return mu;
}


// ���ނ̍ޗ��ƐڐG����������C�������l�������ڐG�͂����߂�D"15_����`�o�l�ƏՓ˂̃��f����"�Q�l�D
// �v�Z����"���{�f���Ȋw�, 2004, 13.4: 233-240."�̔���`�΂˂ւ̎�v�Z�ɂ��K���D 
double Tribology::Tsuji::calc
(
	double k,		// in:	[N/x^1.5]:	����`�����D
	double zeta,	// in:	[-]		 :	������D
	double m,		// in:	[kg]	 :	�������ʁD
	double v,		// in:	[m/s]	 :	�ڋߑ��x�D
	double x		// in:	[m]		 :	�ڋߗʁD
) {
	double c = 2 * zeta * sqrt(1.5 * m * k);
	double Fn = k * pow(x, 1.5) + c * v * pow(x, 0.25);
	// ���l�̂Ƃ���0�ɕ␳�D
	return std::max(Fn, 1e-20);
}

// Kelvin & Voigt�̐ڐG���f����p���Č������l�������ڐG�͂��v�Z 
double Tribology::KelvinVoigt::calc
(
	double k,		// in:	[N/x]:	���`�����D
	double zeta,	// in:	[-]		 :	������D
	double m,		// in:	[kg]	 :	�������ʁD
	double v,		// in:	[m/s]	 :	�ڋߑ��x�D
	double x		// in:	[m]		 :	�ڋߗʁD
) {
	double c = 2 * zeta * sqrt(m * k);
	double Fn = k * x + c * v;
	// ���l�̂Ƃ���0�ɕ␳�D
	return std::max(Fn, 1e-20);
}


// �������BRAIN�Ŏ�������Ă���HamrockDowson�̖��������̎��D
//double Tribology::HamrockDowson_old::calc
//(						// out:	[m]:	���������D
//	double w,			// in:	[N]:	�ڐG�׏d�D
//	double E,			// in:	[Pa]:	���������O���D
//	double Rx,			// in:	[m]:	�]��������ȗ����a�D
//	double Ry,			// in:	[m]:	�������ȗ����a�D
//	double alpha,		// in:	[1/Pa]:	���͔S�x�W���D
//	double eta,			// in:	[Pas]:	���̔S�x�D
//	double u,			// in:	[m/s]:	�]���葬�x�D
//	double lm,			// in:	[m]:	���j�X�J�X�����D�i<-1:�����v�Z�C1~0:���̂܂܊|���Z�C>0:���j�X�J�X�����g�p�D�j
//	double bh			// in:	[m]:	�ڐG�ȉ~�Z�a�D
//) {
//	double W = w / E / Rx / Rx;
//	double G = alpha * E;
//	double U = eta * u / E / Rx;
//
//	double H = 3.63 * pow(U, 0.68) * pow(G, 0.49)
//		* pow(W, -0.073) * (1.0 - exp(-0.7031 * pow(Ry/Rx, 0.636)));
//
//	// Starvation
//	double starv;
//	if (lm <=  -1.0)
//		starv = 1.0;
//
//	else if (lm <=  0.0)
//		starv = fabs(lm);
//
//	else {
//		double f  = 1.0 + lm / bh;
//		double fc = 1.0 + 3.34 * pow(H * Rx * Rx / bh / bh, 0.56);
//		if (f > fc)
//			starv = 1.0;
//		else
//			starv = pow((f - 1.0)/(fc - 1.0), 0.25);
//	}
//	H *= starv;
//
//	double h  =  H * Rx;
//	return h;
//}
//
//// �������BRAIN�Ɏ�������Ă������D�����炭�ԈႢ���Ǝv�����ꉞ�f�ځD
//double Tribology::Aihara_old::calc
//(						// out:	[Nm]:	�g���N�D
//	double Rx,			// in:	[m]:	�����������ȗ��D
//	double Ry,			// in:	[m]:	�]������������ȗ��D
//	double a,			// in:	[m]:	�ڐG�ȉ~���a�D
//	double b,			// in:	[m]:	�ڐG�ȉ~�Z�a�D
//	double w,			// in:	[N]:	�ڐG�׏d�D
//	double E,			// in:	[Pa]:	����0,1�̓��������O���D
//	double um,			// in:	[m/s]:	���ϑ��x�D
//	double eta,			// in:	[Pas]:	���̔S�x�D
//	double alpha,		// in:	[1/Pa]:	���̈��͔S�x�W���D
//	double beta,		// in:	[1/K]:	���̉��x�S�x�W���D
//	double k			// in:	[W/m2K]:���̔M�`�����D
//) {
//	double W = w / E / Ry / Ry;
//	double G = alpha * E;
//	double U = eta * um / E / Ry;
//
//	double L = eta * beta * um * um / k;
//	double p = 3.94 / (3.94 + pow(L, 0.62));
//	double H = (3.63 * pow(U, 0.68) * pow(G, 0.49) * pow(W, -0.073)
//		* (1.0 - exp(-0.7031 * pow(Rx / Ry, 0.636)))) * p;
//
//	double D  = 11.81221 * G * U / W / sqrt(W);
//	double c  = exp((0.658228 + 0.0106383 * pow(L, 0.3174603)) * log(D)
//		- 0.162625 - log(1.0 + 0.411483 * pow(L, 0.707071)));
//	double tr = c * W * Ry / alpha;
//
//	double M = 2.0 * pow(w, 0.3) * tr * Ry * a;
//	return M;
//	//double tr = c * fratio * W * ry / alpha_;
//	// w *= fratio;
//	//double lm			// Meniscus length
//	//double starv;
//	//if (lm <= 0.0)
//	//	starv = 1.0;
//	//else {
//	//	double a1;
//	//	double f = 1.0 + lm / b;
//	//	if ((0.66-0.003*D) > 0.0)
//	//		a1 = 950.0 * pow(D, 0.05) / (1000.0+(0.2*D+9.0)*pow(L, 0.66-0.003*D));
//	//	else
//	//		a1 = 950.0 * pow(D, 0.05) / 1000.0;
//	//	double b  = (pow(1.0+L, 0.06*pow(D, -0.3)-0.035)-0.05*pow(L, 0.0001))*pow(5.0+D, 0.068);
//	//	double f2 = sqrt(f * f - 1.0);
//	//	double s  = (f2 * f - log(f + f2)) / c;
//	//	if (s>0.0) {
//	//		s = 0.6551853*pow(s, 2.0/3.0);
//	//		double s2 = a1*pow(s, b);
//	//		starv = s2 / (1.0 + s2);
//	//	}
//	//	else {
//	//		starv = 0.0;
//	//	}
//	//}
//	//double M = starv * 2.0 * pow(w, 0.3) * tr * ry * a;
//	//if (fratio <= 0.0)
//	//	M = 0.0;
//};
//

