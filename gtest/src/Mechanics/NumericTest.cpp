#include "pch.h"


class NumericTest : public ::testing::Test {
};

TEST_F(NumericTest, Law_of_cosines0) {

	double return_ = Numeric::Law_of_cosines(2.0, 1.0, sqrt(3));
	double answer_ = cos(Numeric::pi / 2);
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
};

TEST_F(NumericTest, Law_of_cosines1) {

	double return_ = Numeric::Law_of_cosines(sqrt(3), 2.0, 1.0);
	double answer_ = cos(Numeric::pi / 3);
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
};

TEST_F(NumericTest, Law_of_cosines2) {

	double return_ = Numeric::Law_of_cosines(sqrt(2), 1.0, 1.0);
	double answer_ = cos(Numeric::pi / 2);
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
};

TEST_F(NumericTest, Law_of_cosines3) {

	double return_ = Numeric::Law_of_cosines(1.0, sqrt(2), 1.0);
	double answer_ = cos(Numeric::pi / 4);
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
};

TEST_F(NumericTest, Square) {
	double return_ = Numeric::Square(3.0);
	double answer_ = 9.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
};

TEST_F(NumericTest, Cube) {
	double return_ = Numeric::Cube(3.0);
	double answer_ = 27.0;
	double delta_ = abs(return_ - answer_);

	EXPECT_GT(1e-2, delta_);
};

TEST_F(NumericTest, Roundoff0) {
	int return0 = Numeric::Roundoff(0.499999999999);
	int answer0 = 0;

	EXPECT_EQ(return0, answer0);

	int return1 = Numeric::Roundoff(0.50000000000);
	int answer1 = 1;

	EXPECT_EQ(return1, answer1);

	int return2 = Numeric::Roundoff(-0.500000000000);
	int answer2 = 0;

	EXPECT_EQ(return2, answer2);

	int return3 = Numeric::Roundoff(-0.500000000001);
	int answer3 = -1;

	EXPECT_EQ(return3, answer3);
};

TEST_F(NumericTest, EffectiveCenter0) {

	double r1 = 10.0;
	double r2 = 15.0;
	double dx = 1.0;

	double return0 = Numeric::EffectiveCenter(r1, r2, dx);

	// éñëOÇÃéËåvéZÇ∆ÇÃî‰ärÅD
	double answer0 = 9.46644333025635;
	EXPECT_EQ(return0, answer0);
};

