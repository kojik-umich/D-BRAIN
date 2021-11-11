/*******************************************************************************
!								"Unit.cpp"
!													2019/11/13	[Core-T]	���
!
!
!	�P�ʕϊ��̂��߂̃I�u�W�F�N�g�D
!
!
!*******************************************************************************/

#include "Unit.h"


// [deg] ���� [rad] �ւ̕ϊ��D
double Unit::deg2rad(double deg) {
	double rad = deg * Numeric::pi / 180;
	return rad;
}

// [rad] ���� [deg] �ւ̕ϊ��D
double Unit::rad2deg(double rad) {
	double deg = rad / Numeric::pi * 180;
	return deg;
}

// [rpm] ���� [rad/s] �ւ̕ϊ��D
double Unit::rpm2rads(double rpm) {
	double rads = rpm * Numeric::pi / 30;
	return rads;
}

// [rpm] ���� [rad/s] �ւ̕ϊ��D
Vector3d Unit::rpm2rads(const Vector3d& rpm) {
	Vector3d rads = rpm * Numeric::pi / 30;
	return rads;
}

// [rad/s] ���� [rpm] �ւ̕ϊ��D
double Unit::rads2rpm(double rads) {
	double rpm = rads / Numeric::pi * 30;
	return rpm;
}

// [rad/s] ���� [rpm] �ւ̕ϊ��D
Vector3d Unit::rads2rpm(const Vector3d&rads) {
	Vector3d rpm = rads / Numeric::pi * 30;
	return rpm;
}

// [s] ���� [ms] �ւ̕ϊ��D
double Unit::s2ms(double s) {
	double ms = 1e3 * s;
	return ms;
}

// [ms] ���� [s] �ւ̕ϊ��D
double Unit::ms2s(double ms) {
	double s = 1e-3 * ms;
	return s;
}

// [m] ���� [mm] �ւ̕ϊ��D
double Unit::m2mm(double m) {
	double mm = 1e3 * m;
	return mm;
}

// [m] ���� [mm] �ւ̕ϊ��D�i�x�N�g���j
Vector3d Unit::m2mm(const Vector3d&m) {
	Vector3d mm = 1e3 * m;
	return mm;
}

// [mm] ���� [m] �ւ̕ϊ��D
double Unit::mm2m(double mm) {
	double m = 1e-3 * mm;
	return m;
}

// [mm] ���� [m] �ւ̕ϊ��D�i�x�N�g���j
Vector3d Unit::mm2m(const Vector3d&mm) {
	Vector3d m = 1e-3 * mm;
	return m;
}

// [Pa*s] ���� [cP] �ւ̕ϊ��D
double Unit::Pas2cP(double Pa) {
	double cP = 1e3 * Pa;
	return cP;
}

// [m/s] ���� [mm/s] �ւ̕ϊ��D
double Unit::ms2mms(double ms) {
	double mms = 1e3 * ms;
	return mms;
}

// [N] ���� [kgf] �ւ̕ϊ��D
double Unit::N2kgf(double N) {
	double kgf = 9.80665 * N;
	return kgf;
}

// [Pa] ���� [kgf/mm2] �ւ̕ϊ��D
double Unit::Pa2kgfmm2(double Pa) {
	double kgfmm2 = 1.0197e-7 * Pa;
	return kgfmm2;
}

// [kgf/mm2] ���� [Pa] �ւ̕ϊ��D
double Unit::kgfmm22Pa(double kgfmm2) {
	double Pa = 9.80665e6 * kgfmm2;
	return Pa;
}

// [1/Pa] ���� [mm2/kgf] �ւ̕ϊ��D
double Unit::Painv2mm2kgf(double Pa) {
	return Unit::kgfmm22Pa(Pa);
}

// [mm2/kgf] ���� [1/Pa] �ւ̕ϊ��D
double Unit::mm2kgf2Painv(double kgfmm2) {
	return Unit::Pa2kgfmm2(kgfmm2);
}

// [GPa] ���� [Pa] �ւ̕ϊ��D
double Unit::GPa2Pa(double GPa) {
	double Pa = 1e9 * GPa;
	return Pa;
}

// [Pa] ���� [GPa] �ւ̕ϊ��D
double Unit::Pa2GPa(double Pa) {
	double GPa = 1e-9 * Pa;
	return GPa;
}

// [Pa] ���� [MPa] �ւ̕ϊ��D
double Unit::Pa2MPa(double Pa) {
	double MPa = 1e-6 * Pa;
	return MPa;
}

// [��m] ���� [m] �ւ̕ϊ��D"��"���łĂȂ��̂�"u"�ő�p�D
double Unit::um2m(double um) {
	double m = 1e-6 * um;
	return m;
}

// [m] ���� [��m] �ւ̕ϊ��D"��"���łĂȂ��̂�"u"�ő�p�D
double Unit::m2um(double m) {
	double um = 1e6 * m;
	return um;
}

// [Nmm] ���� [Nm] �ւ̕ϊ��D
double Unit::Nmm2Nm(double Nmm) {
	double Nm = 1e-3 * Nmm;
	return Nm;
}

// [Nmm] ���� [Nm] �ւ̕ϊ��D
Vector3d Unit::Nmm2Nm(const Vector3d&Nmm) {
	Vector3d Nm = 1e-3 * Nmm;
	return Nm;
}

// [Nm] ���� [Nmm] �ւ̕ϊ��D
double Unit::Nm2Nmm(double Nm) {
	double Nmm = 1e3 * Nm;
	return Nmm;
}

// [Nm] ���� [Nmm] �ւ̕ϊ��D
Vector3d Unit::Nm2Nmm(const Vector3d&Nm) {
	Vector3d Nmm = 1e3 * Nm;
	return Nmm;
}

