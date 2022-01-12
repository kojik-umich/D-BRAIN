/*
この例は、CostFunctionのDynamicAutoDiffCostFunctionバリアントを使用する方法を示しています。
DynamicAutoDiffCostFunctionは、コンパイル時にパラメーターブロックの数またはサイズが不明な
場合に使用することを目的としています。この例では、1次元の廊下を横断するロボットをシミュレートし、
廊下の端のノイズオドメトリ測定値とノイズのある範囲測定値を使用します。この例では、ノイズの多い
オドメトリとセンサーの読み取り値を融合することにより、各タイムステップでのロボットのポーズの
最尤推定（MLE）を計算する方法を示します。ロボットは原点から開始し、「-corridor_length」フラグで
指定された固定長のコリドーの終点まで移動します。一連のモーションコマンドを実行して、
「-pose_separation」フラグで指定された固定長を前方に移動します。どのポーズで、相対的な
走行距離測定値と、廊下の端までの距離の範囲の読み取り値を受け取ります。オドメトリの読み取り値は、
「-odometry_stddev」フラグで指定されたガウスノイズと標準偏差で描画され、範囲の読み取り値は、
「-range-stddev」フラグで指定された標準偏差で同様に描画されます。この問題には2つのタイプの
残差があります。
1）OdometryConstraint残差。これは、ロボットの連続するポーズ推定間のオドメトリ読み取り値を説明します。
2）各ポーズから観測された範囲の読み取り値の誤差を説明するRangeConstraint残差。
OdometryConstraint残余は、固定パラメーターブロックサイズが1の
AutoDiffCostFunctionとしてモデル化されます。これは、解決される相対オドメトリです。
ロボットの連続するポーズのペアの間。観測された相対オドメトリ値と計算された相対オドメトリ値の差は、
オドメトリ測定値の既知の標準偏差によって重み付けされてペナルティが課せられます。
RangeConstraint残余は、DynamicAutoDiffCostFunctionとしてモデル化されます。
これは、相対的な走行距離の推定値を合計して、ロボットの推定されたグローバルポーズを計算し、
次に、予想される範囲の読み取り値を計算します。次に、観測された範囲の読み取り値と予想された
範囲の読み取り値の差が、センサーの読み取り値の標準偏差によって重み付けされてペナルティが
課せられます。ロボットのポーズの数はコンパイル時に不明であるため、このコスト関数は
DynamicAutoDiffCostFunctionとして実装されます。例の出力は、オドメトリと範囲の読み取り値の
初期値、およびロボットのすべてのポーズの範囲とオドメトリのエラーです。MLEを計算した後、
計算されたポーズと修正されたオドメトリ値が、対応する範囲とオドメトリエラーとともに印刷されます。
ノイズの多いシステムのMLEとして、エラーはゼロに減少しませんが、オドメトリの推定値は、
ロボットのすべてのオドメトリと範囲の読み取り値の最尤法を最大化するように更新されることに注意してください。

p_0、..、p_Nを（N + 1）ロボットのポーズとします。ここで、ロボットはp_0から始まり、
p_Nで終わる廊下を移動します。p_0が座標系の原点であると仮定します。オドメトリu_iは、
ポーズp_（i-1）とp_iの間で観測された相対オドメトリであり、範囲読み取り値y_iは、
ポーズp_iからのコリドーの端の範囲読み取り値です。オドメトリと範囲の読み取り値は
どちらもノイズが多いですが、信念が最適化されるように、修正されたオドメトリ値
u * _0からu * _（N-1）の最尤推定（MLE）を計算したいと思います
：Belief（u * _（0：N-1）| u_（0：N-1）、y_（0：N-1））1。
= P（u * _（0：N-1）| u_（0：N-1 ）、y_（0：N-1））2。
\ propto P（y_（0：N-1）| u * _（0：N-1）、u_（0：N-1））P（u * _（0：N-1）| u_（0：N-1））3。
= \ prod_i {P（y_i | u * _（0：i））P（u * _i | u_i）} 4.
ここで、下付き文字「（0：i）」は、その変数のすべてのタイムステップ0からiまでの
エントリを示すための省略形として使用されます。ベイズの定理は、式を導出するために
使用されます。2から3であり、オドメトリ観測と範囲読み取りの独立性は、3から4を
導出するために排除されます。したがって、スケールまでの信念は、各ポーズに2つ、
各ポーズに2つの用語の積として因数分解されます。用語範囲の読み取りには1つの用語
P（y_i | u * _（0：i））と走行距離測定の読み取りには1つの用語P（u * _i | u_i）が
あります。範囲の読み取りの用語はに依存することに注意してください。すべての
オドメトリ値u * _（0：i）、オドメトリ項P（u * _i | u_i）は、単一の値u_iのみに依存します。
範囲読み取り値とオドエムトリー確率項の両方が正規分布としてモデル化されます。
、および形式は次のとおりです。
p（x）\ propto \ exp {-（（x --x_mean）/ x_stddev）^ 2}
ここで、xはMLEオドメトリu *または範囲読み取り値yのいずれかを指し、x_meanは
オドメトリ項の対応する平均値uです。 、およびy_expectedは、以前のすべての
オドメトリ項に基づく予想範囲の読み取り値です。したがって、MLEは、最小化する値
x *を見つけることによって見つけられます。
x* = \ arg \ min {（（x --x_mean）/ x_stddev）^ 2}
は、非線形最小二乗形式であり、Ceresによる解決に適しています。 。
非線形成分は、x_meanの計算から生じます。
Ceresが最適化する残差の残差（（x --x_mean）/ x_stddev）。
前述のように、各ポーズのオドメトリ項は1つの変数のみに依存し、AutoDiffCostFunctionによって
計算されますが、範囲読み取りの項は、以前のすべてのオドメトリ観測に依存します。
#include <math.h>

#include <cstdio>
#include <vector>

#include "ceres/ceres.h"
#include "ceres/dynamic_numeric_diff_cost_function.h"
#include "gflags/gflags.h"
#include "glog/logging.h"

#include <stdlib.h>

#define RAND_MAX 0x7fff


inline double RandDouble() {
	double r = static_cast<double>(rand());
	return r / RAND_MAX;
}

inline double RandNormal() {
	double x1, x2, w;
	do {
		x1 = 2.0 * RandDouble() - 1.0;
		x2 = 2.0 * RandDouble() - 1.0;
		w = x1 * x1 + x2 * x2;
	} while (w >= 1.0 || w == 0.0);

	w = sqrt((-2.0 * log(w)) / w);
	return x1 * w;
}

using ceres::NumericDiffCostFunction;
using ceres::DynamicNumericDiffCostFunction;
using ceres::NumericDiffOptions;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;
using ceres::CENTRAL;
using ceres::TAKE_OWNERSHIP;

using std::min;
using std::vector;

static constexpr double c_len_ = 30.0;
static constexpr double p_sep_ = 0.5;
static constexpr double o_std_ = 0.1;
static constexpr double r_std_ = 0.01;
static constexpr int kStride = 10;

struct BallScrew {

	typedef NumericDiffCostFunction<
		BallScrew,
		CENTRAL,
		1, 1
	> costfunc;

	BallScrew(double odometry_mean, double odometry_stddev)
		: o_mean(odometry_mean), o_std(odometry_stddev) {
	}

	template <typename T>
	bool operator()(const T* const o, T* residual) const {
		*residual = (*o - o_mean) / o_std;
		return true;
	}

	static costfunc* Create(const double odometry_value) {
		return new costfunc(new BallScrew(
			odometry_value, o_std_));
	}

	const double o_mean;
	const double o_std;
};

struct Calculator {
	typedef DynamicNumericDiffCostFunction<Calculator, CENTRAL>
		dyncostfunc;

	Calculator(int pose_index,
		double range_reading,
		double range_stddev,
		double corridor_length)
		: pose_index(pose_index),
		range_reading(range_reading),
		range_stddev(range_stddev),
		corridor_length(corridor_length) {
	}

	template <typename T>
	bool operator()(T const* const* relative_poses, T* residuals) const {
		T global_pose(0);
		for (int i = 0; i <= pose_index; ++i) {
			global_pose += relative_poses[i][0];
		}
		residuals[0] =
			(global_pose + range_reading - corridor_length) / range_stddev;
		return true;
	}

	// セレス問題に追加するのに便利なRangeConstraintからCostFunctionを作成するファクトリーメソッド．
	static dyncostfunc* Create(const int pose_index,
		const double range_reading,
		vector<double>* odometry_values,
		vector<double*>* parameter_blocks) {
		Calculator* constraint =
			new Calculator(pose_index,
				range_reading,
				r_std_,
				c_len_);

		NumericDiffOptions DiffOptions;
		DiffOptions.relative_step_size = 1e-6;
		//DiffOptions.relative_step_size = 1e30;

		dyncostfunc* cost_function = new dyncostfunc(
			constraint,
			TAKE_OWNERSHIP,
			DiffOptions
		);
		// Add all the parameter blocks that affect this constraint.
		parameter_blocks->clear();
		for (int i = 0; i <= pose_index; ++i) {
			parameter_blocks->push_back(&((*odometry_values)[i]));
			cost_function->AddParameterBlock(1);
		}
		cost_function->SetNumResiduals(1);
		return (cost_function);
	}

	const int pose_index;
	const double range_reading;
	const double range_stddev;
	const double corridor_length;
};

void SimulateRobot(
	vector<double>* odometry_values,
	vector<double>* range_readings
) {
	const int num_steps =
		static_cast<int>(ceil(c_len_ / p_sep_));

	// The robot starts out at the origin.
	double robot_location = 0.0;
	for (int i = 0; i < num_steps; ++i) {
		const double actual_odometry_value =
			min(p_sep_, c_len_ - robot_location);
		robot_location += actual_odometry_value;
		const double actual_range =
			c_len_ - robot_location;
		const double o =
			RandNormal() * o_std_ + actual_odometry_value;
		const double r =
			RandNormal() * r_std_ + actual_range;

		odometry_values->push_back(o);
		range_readings->push_back(r);
	}
}

int main(int argc, char** argv) {

	vector<double> o_val;
	vector<double> r_read;
	SimulateRobot(&o_val, &r_read);

	auto print = [](const double& n) { std::cout << " " << n; };

	std::cout << "before:";
	std::for_each(o_val.cbegin(), o_val.cend(), print);
	std::cout << '\n';

	Problem problem;

	for (int i = 0; i < o_val.size(); ++i) {
		// ポーズiからRangeConstraintのDynamicAutoDiffCostFunctionを作成し、追加します。
		vector<double*> param;
		Calculator::dyncostfunc* f =
			Calculator::Create(
				i, r_read[i], &o_val, &param);
		problem.AddResidualBlock(f, NULL, param);

		// ポーズ i の OdometryConstraint 用に AutoDiffCostFunction を作成し、追加します。
		problem.AddResidualBlock(
			BallScrew::Create(o_val[i]),
			NULL,
			&(o_val[i]));
	}

	Solver::Options solver_options;
	solver_options.minimizer_progress_to_stdout = true;

	Solver::Summary summary;
	printf("Solving...\n");
	Solve(solver_options, &problem, &summary);

	printf("Done.\n");
	std::cout << summary.FullReport() << "\n";

	std::cout << "after:";
	std::for_each(o_val.cbegin(), o_val.cend(), print);
	std::cout << '\n';

	return 0;
}

*/

