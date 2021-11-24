#pragma once
#include "Rigid.h"
#include "BS_Cylinder.h"

class BS_Shaft : public Rigid {

public:
	int nCY;			// num-cylinder．現在は１で固定．ねじ軸伸縮考慮で増える？
	BS_Cylinder*CY;		// 螺旋溝円筒部材．
	double tan_thy;		// y軸回りを拘束したときの軸方向とx軸のなす角度θy
	double tan_thz;		// z軸回りを拘束したときの軸方向とx軸のなす角度θz
	// int free_num;		// 拘束されていない自由度の数．[0, 5] の範囲．（軸回転方向は必ず拘束するため最大5）

	struct Memory {
		Vector3d v;
		Vector3d w;
	} mem;

public:
	void Nothing(void);	// 具象クラスであることを示すため，意味のない関数を定義する．
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
