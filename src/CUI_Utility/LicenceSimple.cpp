#include "LicenceSimple.h"

bool LicenceSimple::CheckTime(
	int y0,		// �ݒ莞��
	int m0, 
	int d0, 
	int y1,		// ���ݎ���
	int m1, 
	int d1
) {
	if (y1 > y0)
		return false;

	if (y1 < y0)
		return true;

	if (m1 > m0)
		return false;

	if (m1 < m0)
		return true;

	if (d1 > d0)
		return false;

	return true;
}

bool LicenceSimple::hasTime(int y0, int m0, int d0) {

	time_t now = time(NULL);
	struct tm pnow;
	errno_t error = localtime_s(&pnow, &now);

	int y1 = pnow.tm_year + 1900;
	int m1 = pnow.tm_mon + 1;
	int d1 = pnow.tm_mday;

	bool has = LicenceSimple::CheckTime(y0, m0, d0, y1, m1, d1);

	if (has)
		printf("���C�Z���X�L���������ł��D");
	else
		printf("���C�Z���X�L�������O�ł��D");
	
	printf("%2d�N%2d��%2d���܂ŁD\n", y0, m0, d0);

	return has;
}




