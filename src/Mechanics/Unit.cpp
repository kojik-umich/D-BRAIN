/*******************************************************************************
!								"Unit.cpp"
!													2019/11/13	[Core-T]	楠崎
!
!
!	単位変換のためのオブジェクト．
!
!
!*******************************************************************************/

#include "Unit.h"


// [deg] から [rad] への変換．
double Unit::deg2rad(double deg) {
	double rad = deg * Numeric::pi / 180;
	return rad;
}

// [rad] から [deg] への変換．
double Unit::rad2deg(double rad) {
	double deg = rad / Numeric::pi * 180;
	return deg;
}

// [rpm] から [rad/s] への変換．
double Unit::rpm2rads(double rpm) {
	double rads = rpm * Numeric::pi / 30;
	return rads;
}

// [rpm] から [rad/s] への変換．
Vector3d Unit::rpm2rads(const Vector3d& rpm) {
	Vector3d rads = rpm * Numeric::pi / 30;
	return rads;
}

// [rad/s] から [rpm] への変換．
double Unit::rads2rpm(double rads) {
	double rpm = rads / Numeric::pi * 30;
	return rpm;
}

// [rad/s] から [rpm] への変換．
Vector3d Unit::rads2rpm(const Vector3d&rads) {
	Vector3d rpm = rads / Numeric::pi * 30;
	return rpm;
}

// [s] から [ms] への変換．
double Unit::s2ms(double s) {
	double ms = 1e3 * s;
	return ms;
}

// [ms] から [s] への変換．
double Unit::ms2s(double ms) {
	double s = 1e-3 * ms;
	return s;
}

// [m] から [mm] への変換．
double Unit::m2mm(double m) {
	double mm = 1e3 * m;
	return mm;
}

// [m] から [mm] への変換．（ベクトル）
Vector3d Unit::m2mm(const Vector3d&m) {
	Vector3d mm = 1e3 * m;
	return mm;
}

// [mm] から [m] への変換．
double Unit::mm2m(double mm) {
	double m = 1e-3 * mm;
	return m;
}

// [mm] から [m] への変換．（ベクトル）
Vector3d Unit::mm2m(const Vector3d&mm) {
	Vector3d m = 1e-3 * mm;
	return m;
}

// [Pa*s] から [cP] への変換．
double Unit::Pas2cP(double Pa) {
	double cP = 1e3 * Pa;
	return cP;
}

// [m/s] から [mm/s] への変換．
double Unit::ms2mms(double ms) {
	double mms = 1e3 * ms;
	return mms;
}

// [N] から [kgf] への変換．
double Unit::N2kgf(double N) {
	double kgf = 9.80665 * N;
	return kgf;
}

// [Pa] から [kgf/mm2] への変換．
double Unit::Pa2kgfmm2(double Pa) {
	double kgfmm2 = 1.0197e-7 * Pa;
	return kgfmm2;
}

// [kgf/mm2] から [Pa] への変換．
double Unit::kgfmm22Pa(double kgfmm2) {
	double Pa = 9.80665e6 * kgfmm2;
	return Pa;
}

// [1/Pa] から [mm2/kgf] への変換．
double Unit::Painv2mm2kgf(double Pa) {
	return Unit::kgfmm22Pa(Pa);
}

// [mm2/kgf] から [1/Pa] への変換．
double Unit::mm2kgf2Painv(double kgfmm2) {
	return Unit::Pa2kgfmm2(kgfmm2);
}

// [GPa] から [Pa] への変換．
double Unit::GPa2Pa(double GPa) {
	double Pa = 1e9 * GPa;
	return Pa;
}

// [Pa] から [GPa] への変換．
double Unit::Pa2GPa(double Pa) {
	double GPa = 1e-9 * Pa;
	return GPa;
}

// [Pa] から [MPa] への変換．
double Unit::Pa2MPa(double Pa) {
	double MPa = 1e-6 * Pa;
	return MPa;
}

// [μm] から [m] への変換．"μ"が打てないので"u"で代用．
double Unit::um2m(double um) {
	double m = 1e-6 * um;
	return m;
}

// [m] から [μm] への変換．"μ"が打てないので"u"で代用．
double Unit::m2um(double m) {
	double um = 1e6 * m;
	return um;
}

// [Nmm] から [Nm] への変換．
double Unit::Nmm2Nm(double Nmm) {
	double Nm = 1e-3 * Nmm;
	return Nm;
}

// [Nmm] から [Nm] への変換．
Vector3d Unit::Nmm2Nm(const Vector3d&Nmm) {
	Vector3d Nm = 1e-3 * Nmm;
	return Nm;
}

// [Nm] から [Nmm] への変換．
double Unit::Nm2Nmm(double Nm) {
	double Nmm = 1e3 * Nm;
	return Nmm;
}

// [Nm] から [Nmm] への変換．
Vector3d Unit::Nm2Nmm(const Vector3d&Nm) {
	Vector3d Nmm = 1e3 * Nm;
	return Nmm;
}

