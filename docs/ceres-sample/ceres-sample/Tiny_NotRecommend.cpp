/*
#include "ceres/ceres.h"
#include "ceres/tiny_solver.h"


// これから使う予定
#include "ceres/numeric_diff_cost_function.h"
using ceres::CENTRAL;
using ceres::NumericDiffCostFunction;
using ceres::Problem;



using Eigen::VectorXd;

class ExampleAllDynamic {
public:
	typedef double Scalar;
	enum {
		NUM_RESIDUALS = Eigen::Dynamic,
		NUM_PARAMETERS = Eigen::Dynamic,
	};

	int NumResiduals() const {
		return 2;
	}

	int NumParameters() const {
		return 3;
	}

	bool operator()(const Scalar* parameters,
		Scalar* residuals,
		Scalar* jacobian) const
	{
		Scalar x = parameters[0];
		Scalar y = parameters[1];
		Scalar z = parameters[2];

		residuals[0] = x + 2 * y + 4 * z;
		residuals[1] = y * z;

		if (jacobian) {
			std::cout << "ヤコビアンが計算されています" << std::endl << std::endl;

			jacobian[0 * 2 + 0] = 1;
			jacobian[0 * 2 + 1] = 0;

			jacobian[1 * 2 + 0] = 2;
			jacobian[1 * 2 + 1] = z;

			jacobian[2 * 2 + 0] = 4;
			jacobian[2 * 2 + 1] = y;
		}
		return true;
	}
};

void check_variables(VectorXd&x, VectorXd&residuals, ExampleAllDynamic&f) {

	f(x.data(), residuals.data(), NULL);
	std::cout << "model function f(x) = " << residuals.squaredNorm() / 2.0 << std::endl;
	std::cout << "x0 = " << std::endl << x << std::endl << std::endl;
	return;
}

int main(int argc, char** argv) {

	VectorXd x(3);
	x << 0.76026643, -30.01799744, 0.55192142;

	ExampleAllDynamic f;

	VectorXd residuals(2);

	check_variables(x, residuals, f);

	ceres::TinySolver<ExampleAllDynamic> solver;
	solver.Solve(f, &x);

	check_variables(x, residuals, f);

	return 0;
}


*/
