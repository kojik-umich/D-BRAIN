#include <cmath>
#include <iostream>
#include <vector>
#include <ceres/ceres.h>
#include "BS_BallScrew.h"

struct BS_CostFunctor {

private:
	std::shared_ptr<BS_BallScrew> BS_;
	double v0_;
	double w0_;

public:
	explicit BS_CostFunctor(
		std::shared_ptr<BS_BallScrew> BS,
		double v0,
		double w0
	) : BS_(BS), v0_(v0), w0_(w0) {
	}
	virtual ~BS_CostFunctor() = default;

	bool operator()(double const* const* parameters, double* residual) const {

		this->BS_->set_y0(parameters[0], this->v0_, this->w0_);
		this->BS_->get_F1(residual);

		return true;
	}

	static auto Create(
		std::shared_ptr<BS_BallScrew> BS,
		double v0,
		double w0,
		int n,
		const ceres::NumericDiffOptions&diffoptions
	) {
		auto cost_function =
			new ceres::DynamicNumericDiffCostFunction<
			BS_CostFunctor,
			ceres::CENTRAL
			>(
				new BS_CostFunctor(BS, v0, w0),
				ceres::TAKE_OWNERSHIP,
				diffoptions);

		cost_function->AddParameterBlock(n);
		cost_function->SetNumResiduals(n);

		return cost_function;
	}
};
