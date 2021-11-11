#pragma once
#include "Rigid.h"
#include "BS_Cylinder.h"

class BS_Shaft : public Rigid {

public:
	int nCY;			// num-cylinderDŒ»İ‚Í‚P‚ÅŒÅ’èD‚Ë‚¶²Lkl—¶‚Å‘‚¦‚éH
	BS_Cylinder*CY;		// —†ùa‰~“›•”ŞD
	double tan_thy;		// y²‰ñ‚è‚ğS‘©‚µ‚½‚Æ‚«‚Ì²•ûŒü‚Æx²‚Ì‚È‚·Šp“xƒÆy
	double tan_thz;		// z²‰ñ‚è‚ğS‘©‚µ‚½‚Æ‚«‚Ì²•ûŒü‚Æx²‚Ì‚È‚·Šp“xƒÆz

public:
	void Nothing(void);	// ‹ïÛƒNƒ‰ƒX‚Å‚ ‚é‚±‚Æ‚ğ¦‚·‚½‚ßCˆÓ–¡‚Ì‚È‚¢ŠÖ”‚ğ’è‹`‚·‚éD
	void allocate(const std::vector<BS_In::Cylinder>& cylinders);
	void init(const std::vector<BS_In::Cylinder>& cylinders, const bool(&v_const)[3], const bool(&w_const)[3], double tan_thy, double tan_thz, double v0, double w0);
	void init_pos(double v0, double w0);
	void get_y0(double * y0);
	void set_y0(const double * y0, double v0, double w0);
	void set_dx(void);
	void set_y_(const Vector3d & x, const Vector3d & v, const Quaterniond & q, const Vector3d & w);
	void get_dydt_(const Vector3d & F, const Vector3d & T, double dvdt, double dwdt, double * dydt);
	Vector4d get_dqdt_(bool wy_const, bool wz_const, double tan_thy, double tan_thz);
	Vector3d get_dwdt_(const Vector3d&T, bool wy_const, bool wz_const, double dwdt0);
	BS_Shaft();
	~BS_Shaft();
};
