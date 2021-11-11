#include "pch.h"



// トラクション係数の数値として常に0.2を出力．
double Tribology::Stab_Traction::calc
(						// out:	[-]:	トラクション係数．
	double eta0,		// in:	[Pas]:	油の粘度．
	double p0,			// in:	[Pa]:	最大面圧．
	double v0,			// in:	[m/s]:	滑り速度．
	double u0			// in:	[m/s]:	転がり速度．
) {
	return 0.2;
}

// クーロン摩擦係数の数値として常に0.3を出力．ただし，滑り速度が0に近い場合，0を出力．
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