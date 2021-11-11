#pragma once
#include <math.h>

namespace Numeric {
	const double pi= 3.14159265358979323846;				// ‰~Žü—¦ƒÎ
	const double ot= 0.33333333333333333333;				// 3•ª‚Ì1
	const double pi_half= 1.57079632679489661923;			// ƒÎ / 2
	const double pi_quarter_inv= 1.2732395447351628;		// 4 / ƒÎ

	double Law_of_cosines(double a, double b, double c);
	double Square(double x);
	double Cube(double x);
	int Roundoff(double x);
	double EffectiveCenter(double r1, double r2, double dx);
};

