#pragma once
#include "pch.h"

class bal_StabCage : public bal_Cage {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	void init(const B4P_In&FI) {
		return;
	};
	// �ڐG�_���Ƃ���2�_�o�́i1�_�ڂ��ʂ��\�������C2�_�ڂ��ʂɐH�����ށj
	int get_ContactPoint(const Vector3d & BL_x, double BL_r, int np, Vector3d*x, double*k, int*ptt) {
		x[1] = Vector3d(0, 5.0e-3, 30.0e-3);
		k[1] = 1e5;
		x[0] = Vector3d(0, 0, 30.0);
		k[0] = 1e5;
		return 2;
	};
};