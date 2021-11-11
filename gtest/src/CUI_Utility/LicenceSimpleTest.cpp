#include "pch.h"

class LicenceSimpleTest : public ::testing::Test {
};

TEST_F(LicenceSimpleTest, CheckTime0) {

	// ���̐ݒ莞���D
	int y0 = 2020;
	int m0 = 7;
	int d0 = 15;

	EXPECT_EQ(true, LicenceSimple::CheckTime(y0, m0, d0, 2019, 12, 31));
	EXPECT_EQ(false, LicenceSimple::CheckTime(y0, m0, d0, 2021, 12, 31));
	EXPECT_EQ(true, LicenceSimple::CheckTime(y0, m0, d0, 2020, 6, 31));
	EXPECT_EQ(false, LicenceSimple::CheckTime(y0, m0, d0, 2020, 8, 31));
	EXPECT_EQ(true, LicenceSimple::CheckTime(y0, m0, d0, 2020, 7, 14));
	EXPECT_EQ(false, LicenceSimple::CheckTime(y0, m0, d0, 2020, 7, 16));
	
	return;
};


TEST_F(LicenceSimpleTest, hasTime0) {

	// �͂邩�ߋ��̐ݒ�Ŏ��s���邩�̃e�X�g�D
	EXPECT_EQ(false, LicenceSimple::hasTime(1950, 12, 31));

	// �͂邩�����̐ݒ�Ő������邩�̃e�X�g�D
	EXPECT_EQ(true, LicenceSimple::hasTime(3000, 7, 16));

	return;
};
