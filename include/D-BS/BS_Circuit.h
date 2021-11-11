#pragma once
#include "Ball.h"
#include "BS_In.h"
#include "BS_Cylinder.h"
#include <vector>
using Eigen::Matrix3d;

class BS_Circuit {

public:
	int nBL;		// 回路に含まれている玉数．
	Ball*BL;		// 回路に含まれている玉．
	BS_Cylinder*CY;	// 所属している円筒．おそらくナット．
	int iSP;		// 所属する条番号．
	double th0;		// 回路の-x側位相角．
	double th1;		// 回路の+x側位相角．

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

