#include "Numeric.h"


// �]���藝�̌����D
double Numeric::Law_of_cosines(double a, double b, double c) {
	double cos_th = (b*b + c * c - a * a) / (2 * b*c);
	return cos_th;
}

// 2�悷��D
double Numeric::Square(double x) {
	return x * x;
}

// x��3��
double Numeric::Cube(double x) {
	return x * x * x;
}

// �������ʂ��l�̌ܓ����Cint�^�ɂ��ďo��
int Numeric::Roundoff(double x) {
	int i = (int)(x + 0.5);
	if (x < -0.5)
		return i - 1;
	return i;
}

// 2�̉~�ʂ̓��ߋȗ����a�̐����C2�~�ʂ̌�_��ʂ�悤�Ȍʂ̒��S�ʒu�����߂郁�\�b�h�D
double Numeric::EffectiveCenter
(					// out:	[m]:	r1 �̒��S���猩���C�����~�ʂ̐�[�܂ł̋����D
	double r1,		// in:	[m]:	��ƂȂ�~�ʁD
	double r2,		// in:	[m]:	������̉~�ʁDr1 < r2 �̊֌W�D
	double dx		// in:	[m]:	2�~�ʂ̐ڋߗʂ̍ő�l�Dr1���S��r2���S�����񂾒�����ő���D
	) 
{
	double Rm = 2 * (r1 * r2) / (r1 + r2);	// �ڐG�ʂ̋ȗ����a
	double O2 = r2 + dx - r1;
	double y1 = (Square(r2) - Square(O2) - Square(r1)) / (2 * O2);
	double x1sq = Square(r1) - Square(y1);
	double x0 = y1 + Rm - sqrt(Square(Rm) - x1sq);

	return x0;
}