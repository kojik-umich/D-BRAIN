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
double Tribology::Stab_RollingResist::calc
(						// out:	[Nm]:	トルク．
	double Rx,			// in:	[m]:	転がり方向等価曲率．
	double a,			// in:	[m]:	接触楕円長径．
	double b,			// in:	[m]:	接触楕円短径．
	double w,			// in:	[N]:	接触荷重．
	double E,			// in:	[Pa]:	物体0,1の等価ヤング率．
	double um,			// in:	[m/s]:	平均速度．
	double eta,			// in:	[Pas]:	油の粘度．
	double alpha,		// in:	[1/Pa]:	油の圧力粘度係数．
	double beta,		// in:	[1/K]:	油の温度粘度係数．
	double k			// in:	[W/m2K]:油の熱伝導率．
) {
	return 0.1;
}