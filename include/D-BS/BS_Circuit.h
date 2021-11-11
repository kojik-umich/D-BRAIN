#pragma once
#include "Ball.h"
#include "BS_In.h"
#include "BS_Cylinder.h"
#include <vector>
using Eigen::Matrix3d;

class BS_Circuit {

public:
	int nBL;		// ��H�Ɋ܂܂�Ă���ʐ��D
	Ball*BL;		// ��H�Ɋ܂܂�Ă���ʁD
	BS_Cylinder*CY;	// �������Ă���~���D�����炭�i�b�g�D
	int iSP;		// ���������ԍ��D
	double th0;		// ��H��-x���ʑ��p�D
	double th1;		// ��H��+x���ʑ��p�D

public:
	BS_Circuit(void);
	~BS_Circuit(void);
	void allocate(int n);
	void link(BS_Cylinder*CY, int is);
	void init(const BS_In::Circuit&circuit);
	void get_x0s(Vector3d * xs);
	Vector3d get_t(int ib);
	double get_nd(void);
	double get_r(void);
};

