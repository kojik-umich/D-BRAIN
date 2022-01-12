/*

// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2021 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: sameeragarwal@google.com (Sameer Agarwal)

Levenberg-Marquardtを使用した小さな最小二乗ソルバーで、低レイテンシと低オーバーヘッドの小さな密度の高い問題を解決することを目的としています。実装では、すべての割り当てを前もって行うように注意しているため、解決中にメモリが割り当てられることはありません。これは、多くの同様の問題を解決するときに特に役立ちます。たとえば、グリッド上のすべてのピクセルの逆ピクセル歪み。注：このコードには、Ceresの他の部分を含め、Eigen以外の依存関係はありません。したがって、このファイルを単独で取得して、Ceresの残りの部分なしで別のプロジェクトに配置することができます。以下に基づくアルゴリズム：[1] K. Madsen、H。Nielsen、O。Tingleoff。非線形最小二乗問題の方法。

#include "ceres/ceres.h"
#include "glog/logging.h"
#include "ceres/types.h"

using ceres::CENTRAL;
using ceres::CostFunction;
using ceres::NumericDiffCostFunction;
using ceres::DynamicNumericDiffCostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;
using Eigen::VectorXd;
using std::vector;

class MyCostFunctor {
public:
	bool operator()(double const* const* parameters, double* residuals) const {
		const double* params0 = parameters[0];
		int r = 0;
		for (int i = 0; i < 10; ++i) {
			residuals[r++] = i - params0[i];
			residuals[r++] = params0[i] - i;
		}

		double c_residual = 0.0;
		for (int i = 0; i < 10; ++i) {
			c_residual += pow(params0[i], 2) - 8.0 * params0[i];
		}

		const double* params1 = parameters[1];
		for (int i = 0; i < 5; ++i) {
			c_residual += params1[i];
		}
		residuals[r++] = c_residual;
		return true;
	}
};

int main(int argc, char** argv) {

	int num_residuals = 21;

	vector<double> param_block_0(10, 0.0);
	for (int i = 0; i < 10; ++i)
		param_block_0[i] = 2 * i;

	vector<double> param_block_1(5, 0.0);

	DynamicNumericDiffCostFunction<MyCostFunctor> cost_function(
		new MyCostFunctor());
	cost_function.AddParameterBlock(param_block_0.size());
	cost_function.AddParameterBlock(param_block_1.size());
	cost_function.SetNumResiduals(21);

	ceres::Problem problem;

	ceres::NumericDiffOptions options;
	options.relative_step_size = 1e-6;

	CostFunction* cost_function
		= new NumericDiffCostFunction<
		MyCostFunctor,
		CENTRAL,
		ceres::DYNAMIC,
		num_residuals
		>(
			new MyCostFunctor(),
			ceres::TAKE_OWNERSHIP,
			num_residuals,
			options
			);


	// Run the solver!
	Solver::Options SolverOptions;

	SolverOptions.minimizer_progress_to_stdout = true;
	Solver::Summary summary;

	Solve(SolverOptions, &problem, &summary);

	return 0;
}

//VectorXd x(3);
//x << 0.76026643, -30.01799744, 0.55192142;

//const VectorXd initial_x = x;

//// Build the problem.
//Problem problem;

//ceres::NumericDiffOptions DiffOptions;
//DiffOptions.relative_step_size = 1e-6;

//// Set up the only cost function (also known as residual). This uses
//// numeric differentiation to obtain the derivative (jacobian).
//int num_residuals = 2;
//int num_parameters = 3;

//CostFunction* cost_function =
//	new ceres::DynamicNumericDiffCostFunction<ExampleAllDynamic //,
//	//CENTRAL,
//	//2,	// Dimensions of residuals
//	//3 	// Dimensions of parameters
//	>(
//		new ExampleAllDynamic,
//		ceres::TAKE_OWNERSHIP,
//		num_residuals,
//		DiffOptions
//		);

//problem.AddResidualBlock(cost_function, NULL, x.data());
//problem.SetNumResiduals(21);

//// Run the solver!
//Solver::Options SolverOptions;

//SolverOptions.minimizer_progress_to_stdout = true;
//Solver::Summary summary;

////Solve(SolverOptions, &problem, &summary);

//ExampleAllDynamic f;
//ceres::TinySolver<ExampleAllDynamic> solver;
//solver.Solve(f, &x);

//std::cout << summary.BriefReport() << "\n";
//std::cout << "x : \n" << initial_x << "\n -> \n" << x << "\n";




//ceres::NumericDiffOptions options;
//options.relative_step_size = 1e-6;


*/