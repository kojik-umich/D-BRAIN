#include <cmath>
#include <iostream>
#include <vector>
#include <ceres/ceres.h>

class SpringSystem {

private:
	static double Hertz(double k, double dx) {
		return k * (
			(dx > 0)
			? std::pow(dx, 1.5)
			: std::pow(-dx, 1.5) * 0e-4
			);
	};
	const double k0 = 1.0;
	const double k1 = 2.0;
	const double k2 = 4.0;
	const double F0 = 1.0;

	double x0;
	double x1;

public:
	SpringSystem() {
		this->x0 = 0;
		this->x1 = 0;
	};

	virtual ~SpringSystem() = default;
	SpringSystem& operator=(SpringSystem&& o) = default;

	void getPosition(std::vector<double>& x) const {
		x[0] = this->x0;
		x[1] = this->x1;
	};

	void setPosition(double const* const x) {
		this->x0 = x[0];
		this->x1 = x[1];
	};

	void getForce(double * F) {
		double dx0 = this->x0;
		double dx1 = this->x0 - this->x1;
		double dx2 = this->x1;

		F[0] = -Hertz(this->k0, dx0) - Hertz(this->k1, dx1) + this->F0;
		F[1] = Hertz(this->k1, dx1) - Hertz(this->k2, dx2);
	};
};

struct SpringSystemCost : public ceres::CostFunction {

private:
	std::shared_ptr<SpringSystem> SS_;

public:
	explicit SpringSystemCost(std::shared_ptr<SpringSystem> SS)
		: SS_(SS) {
	}

	bool Evaluate(
		double const* const* parameters,
		double* residuals,
		double** jacobians) const {

		SS_->setPosition(parameters[0]);
		SS_->getForce(residuals);

		return true;
	}

	static auto Create(
		std::shared_ptr<SpringSystem> SS,
		const ceres::NumericDiffOptions&diffoptions
	) {
		auto cost_function =
			new ceres::DynamicNumericDiffCostFunction<
			SpringSystemCost,
			ceres::CENTRAL
			>(
				new SpringSystemCost(SS),
				ceres::TAKE_OWNERSHIP,
				diffoptions);

		int abs_point_size = 2 * 1;
		int residual_size = 7 - 5;

		cost_function->AddParameterBlock(abs_point_size);
		cost_function->SetNumResiduals(residual_size);

		return cost_function;
	}
};

void printXF(std::shared_ptr<SpringSystem> SS);

int main(int argc, char** argv) {

	auto SS = std::make_shared<SpringSystem>();

	const int n = 1 + 1;
	std::vector<double> x(n);
	SS->getPosition(x);
	printXF(SS);

	ceres::NumericDiffOptions diffoptions;
	diffoptions.relative_step_size = 1e-9;

	auto cost = SpringSystemCost::Create(SS, diffoptions);

	ceres::Problem problem;
	problem.AddResidualBlock(cost, NULL, x.data());

	ceres::Solver::Options options;
	options.max_num_iterations = 100;
	options.linear_solver_type = ceres::DENSE_QR;
	options.minimizer_progress_to_stdout = true;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);

	std::cout << summary.FullReport() << "\n";

	printXF(SS);

	return 0;
}

void printXF(std::shared_ptr<SpringSystem> SS) {
	std::vector<double> x(2);
	SS->getPosition(x);
	std::cout << "x =\t" << x[0] << ",\t" << x[1] << "\n";

	SS->getForce(x.data());
	std::cout << "F =\t" << x[0] << ",\t" << x[1] << "\n\n";
}

















//std::cout << "x =\t" << x[0] << ",\t" << x[1] << "\n";

//double f[2];
//SS->getForce(f);
//std::cout << "F =\t" << f[0] << ",\t" << f[1] << "\n";
