#pragma once
#include "Rigid.h"
#include "BS_Cylinder.h"

class BS_Shaft : public Rigid {

public:
	int nCY;			// num-cylinder�D���݂͂P�ŌŒ�D�˂����L�k�l���ő�����H
	BS_Cylinder*CY;		// �����a�~�����ށD
	double tan_thy;		// y�������S�������Ƃ��̎�������x���̂Ȃ��p�x��y
	double tan_thz;		// z�������S�������Ƃ��̎�������x���̂Ȃ��p�x��z
	// int free_num;		// �S������Ă��Ȃ����R�x�̐��D[0, 5] �͈̔́D�i����]�����͕K���S�����邽�ߍő�5�j

	struct Memory {
		Vector3d v;
		Vector3d w;
	} mem;

public:
	void Nothing(void);	// ��ۃN���X�ł��邱�Ƃ��������߁C�Ӗ��̂Ȃ��֐����`����D
	void allocate(const std::vector<BS_In::Cylinder>& cylinders);
	void init(const std::vector<BS_In::Cylinder>& cylinders, const bool(&v_const)[3], const bool(&w_const)[3], double tan_thy, double tan_thz, double v0, double w0);
	void init_pos(double v0, double w0);
	void get_y0(double * y0);
	void set_y0(const double * y0, double v0, double w0);
	void set_dx(void);
	void set_y_(const Vector3d & x, const Vector3d & v, const Quaterniond & q, const Vector3d & w);
	void init_dyn0(void);
	void deinit_dyn0(void);
	void set_dyn_y0(const double * y0);
	void get_dyn_y0(double * y0);
	void get_dydt_(const Vector3d & F, const Vector3d & T, double dvdt, double dwdt, double * dydt);
	Vector4d get_dqdt_(bool wy_const, bool wz_const, double tan_thy, double tan_thz);
	Vector3d get_dwdt_(const Vector3d&T, bool wy_const, bool wz_const, double dwdt0);
	BS_Shaft();
	~BS_Shaft();
};
