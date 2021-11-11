#pragma once
#include <math.h>
#include "Numeric.h"
#include <Eigen\Dense>
using Eigen::Vector3d;

namespace Unit {

	double deg2rad(double deg);
	double rad2deg(double rad);
	double rpm2rads(double rpm);
	Vector3d rpm2rads(const Vector3d & rpm);
	double rads2rpm(double rads);
	Vector3d rads2rpm(const Vector3d & rads);
	double s2ms(double s);
	double ms2s(double ms);
	double m2mm(double m);
	double mm2m(double mm);
	Vector3d m2mm(const Vector3d&m);
	Vector3d mm2m(const Vector3d&mm);
	double Pas2cP(double Pa);
	double ms2mms(double ms);
	double N2kgf(double N);
	double Pa2kgfmm2(double Pa);
	double kgfmm22Pa(double kgfmm2);
	double Painv2mm2kgf(double Pa);
	double mm2kgf2Painv(double kgfmm2);
	double GPa2Pa(double GPa);
	double Pa2GPa(double Pa);
	double Pa2MPa(double Pa);
	double um2m(double um);
	double m2um(double m);
	double Nmm2Nm(double Nmm);
	Vector3d Nmm2Nm(const Vector3d & Nmm);
	double Nm2Nmm(double Nm);
	Vector3d Nm2Nmm(const Vector3d & Nm);
};

