#include "Numeric.h"


// 余弦定理の公式．
double Numeric::Law_of_cosines(double a, double b, double c) {
	double cos_th = (b*b + c * c - a * a) / (2 * b*c);
	return cos_th;
}

// 2乗する．
double Numeric::Square(double x) {
	return x * x;
}

// xの3乗
double Numeric::Cube(double x) {
	return x * x * x;
}

// 小数第一位を四捨五入し，int型にして出力
int Numeric::Roundoff(double x) {
	int i = (int)(x + 0.5);
	if (x < -0.5)
		return i - 1;
	return i;
}

// 2つの円弧の透過曲率半径の線が，2円弧の交点を通るような弧の中心位置を求めるメソッド．
double Numeric::EffectiveCenter
(					// out:	[m]:	r1 の中心から見た，等価円弧の先端までの距離．
	double r1,		// in:	[m]:	基準となる円弧．
	double r2,		// in:	[m]:	もう一つの円弧．r1 < r2 の関係．
	double dx		// in:	[m]:	2円弧の接近量の最大値．r1中心とr2中心を結んだ直線上で測定．
	) 
{
	double Rm = 2 * (r1 * r2) / (r1 + r2);	// 接触面の曲率半径
	double O2 = r2 + dx - r1;
	double y1 = (Square(r2) - Square(O2) - Square(r1)) / (2 * O2);
	double x1sq = Square(r1) - Square(y1);
	double x0 = y1 + Rm - sqrt(Square(Rm) - x1sq);

	return x0;
}