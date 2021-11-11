#include "pch.h"


class UnitTest : public ::testing::Test {
};

TEST_F(UnitTest, deg2rad) {

	double return_ = Unit::deg2rad(150.0);
	double answer_ = Numeric::pi *5 / 6;
	double delta_  = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, rad2deg) {

	double return_ = Unit::rad2deg(Numeric::pi * 5 / 6);
	double answer_ = 150.0;
	double delta_  = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, rpm2rads_d) {

	double return_ = Unit::rpm2rads(30.0);
	double answer_ = Numeric::pi;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, rpm2rads_v) {
	Vector3d input = Vector3d(60, 30, -60);
	Vector3d return_ = Unit::rpm2rads(input);
	Vector3d answer_ = Numeric::pi * Vector3d(2, 1, -2);
	double delta_ = (return_ - answer_).norm();

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, rads2rpm_d) {

	double return_ = Unit::rads2rpm(Numeric::pi);
	double answer_ = 30.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, rads2rpm_v) {

	Vector3d return_ = Unit::rads2rpm(Vector3d(Numeric::pi, Numeric::pi, Numeric::pi));
	Vector3d answer_ = Vector3d(30.0, 30.0, 30.0);
	double delta_ = (return_ - answer_).norm();

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, s2ms) {

	double return_ = Unit::s2ms(1.0);
	double answer_ = 1000.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, ms2s) {

	double return_ = Unit::ms2s(1000.0);
	double answer_ = 1.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, m2mm_d) {

	double return_ = Unit::m2mm(1.0);
	double answer_ = 1000.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, m2mm_v) {

	Vector3d m_      = Vector3d(1.0, 2.0, 3.0);
	Vector3d return_ = Unit::m2mm(m_);
	Vector3d answer_ = Vector3d(1000.0, 2000.0, 3000.0);
	double delta_    = (return_ - answer_).norm();

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, mm2m_d) {

	double return_ = Unit::mm2m(1000.0);
	double answer_ = 1.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, mm2m_v) {

	Vector3d m_ = Vector3d(1000.0, 2000.0, 3000.0);
	Vector3d return_ = Unit::mm2m(m_);
	Vector3d answer_ = Vector3d(1.0, 2.0, 3.0);
	double delta_ = (return_ - answer_).norm();

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, Pas2cP) {

	double return_ = Unit::Pas2cP(1e-3);
	double answer_ = 1.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, ms2mms) {

	double return_ = Unit::ms2mms(1.0);
	double answer_ = 1000.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, N2kgf) {

	double return_ = Unit::N2kgf(5.0);
	double answer_ = 49.0332;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, Pa2kgfmm2) {

	double return_ = Unit::Pa2kgfmm2(3000.0);
	double answer_ = 0.00030591;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, kgfmm22Pa) {

	double return_ = Unit::kgfmm22Pa(0.10);
	double answer_ = 980665;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, Painv2mm2kgf) {

	double return_ = Unit::Painv2mm2kgf(0.10);
	double answer_ = 980665;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, mm2kgf2Painv) {

	double return_ = Unit::mm2kgf2Painv(980665);
	double answer_ = 0.10;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, GPa2Pa) {

	double return_ = Unit::GPa2Pa(1.0);
	double answer_ = 1.0e9;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, Pa2GPa) {

	double return_ = Unit::Pa2GPa(1e9);
	double answer_ = 1.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, Pa2MPa) {

	double return_ = Unit::Pa2MPa(1e6);
	double answer_ = 1.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, um2m) {

	double return_ = Unit::um2m(1);
	double answer_ = 1e-6;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, m2um) {

	double return_ = Unit::m2um(1);
	double answer_ = 1e6;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, Nmm2Nm_d) {

	double return_ = Unit::Nmm2Nm(1);
	double answer_ = 1e-3;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, Nmm2Nm_v) {

	Vector3d return_ = Unit::Nmm2Nm(Vector3d(1.0, 1.0, 1.0));
	Vector3d answer_ = Vector3d(1e-3, 1e-3, 1e-3);
	double delta_ = (return_ - answer_).norm();

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, Nm2Nmm_d) {

	double return_ = Unit::Nm2Nmm(1);
	double answer_ = 1e3;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
	return;
};

TEST_F(UnitTest, Nm2Nmm_v) {

	Vector3d return_ = Unit::Nm2Nmm(Vector3d(1.0, 1.0, 1.0));
	Vector3d answer_ = Vector3d(1e3, 1e3, 1e3);
	double delta_ = (return_ - answer_).norm();

	EXPECT_GT(1e-2, delta_);
	return;
};





























