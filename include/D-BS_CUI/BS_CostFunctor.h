#include <cmath>
#include <iostream>
#include <vector>
#include <ceres/ceres.h>
#include "BS_BallScrew.h"

struct BS_CostFunctor {

private:
	std::shared_ptr<BS_BallScrew> BS_;

public:
	explicit BS_CostFunctor(std::shared_ptr<BS_BallScrew> BS) : BS_(BS) {
	}
	virtual ~BS_CostFunctor() = default;

	bool operator()(double const* const* parameters, double* residual) const {
		return true;
	}

	static auto Create(
		std::shared_ptr<BS_BallScrew> BS,
		const ceres::NumericDiffOptions&diffoptions
	) {
		auto cost_function =
			new ceres::DynamicNumericDiffCostFunction<
			BS_CostFunctor,
			ceres::CENTRAL
			>(
				new BS_CostFunctor(BS),
				ceres::TAKE_OWNERSHIP,
				diffoptions);

		int abs_point_size = 2 * 1;
		int residual_size = 7 - 5;

		cost_function->AddParameterBlock(abs_point_size);
		cost_function->SetNumResiduals(residual_size);

		return cost_function;
	}
};
