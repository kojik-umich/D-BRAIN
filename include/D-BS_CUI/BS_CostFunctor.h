#include <cmath>
#include <iostream>
#include <vector>
#include <ceres/ceres.h>
#include "BS_BallScrew.h"

struct BS_CostSimpleCoulomb {

private:
	std::shared_ptr<BS_BallScrew> BS_;
	double v0_;
	double w0_;

public:

	explicit BS_CostSimpleCoulomb(
		std::shared_ptr<BS_BallScrew> BS,
		double v0,
		double w0
		) : BS_(BS), v0_(v0), w0_(w0) {
	}
	virtual ~BS_CostSimpleCoulomb() = default;

	bool operator()(
		double const* const* parameters,
		double* residual
		) const {

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
			BS_CostSimpleCoulomb,
			ceres::CENTRAL
			>(
				new BS_CostSimpleCoulomb(BS, v0, w0),
				ceres::TAKE_OWNERSHIP,
				diffoptions);

		cost_function->AddParameterBlock(n);
		cost_function->SetNumResiduals(n);

		return cost_function;
	}
};

struct BS_CostFullParameter {

private:
	std::shared_ptr<BS_BallScrew> BS_;
	double v0_;
	double w0_;

public:

	explicit BS_CostFullParameter(
		std::shared_ptr<BS_BallScrew> BS,
		double v0,
		double w0
	) : BS_(BS), v0_(v0), w0_(w0) {
	}
	virtual ~BS_CostFullParameter() = default;

	bool operator()(
		double const* const* parameters,
		double* residual
		) const {

		this->BS_->set_y2(parameters[0], this->v0_, this->w0_);
		this->BS_->get_F2(residual);
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
			BS_CostFullParameter,
			ceres::CENTRAL
			>(
				new BS_CostFullParameter(BS, v0, w0),
				ceres::TAKE_OWNERSHIP,
				diffoptions);

		cost_function->AddParameterBlock(n);
		cost_function->SetNumResiduals(n);

		return cost_function;
	}
};
