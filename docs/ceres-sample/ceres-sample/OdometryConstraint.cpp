/*
���̗�́ACostFunction��DynamicAutoDiffCostFunction�o���A���g���g�p������@�������Ă��܂��B
DynamicAutoDiffCostFunction�́A�R���p�C�����Ƀp�����[�^�[�u���b�N�̐��܂��̓T�C�Y���s����
�ꍇ�Ɏg�p���邱�Ƃ�ړI�Ƃ��Ă��܂��B���̗�ł́A1�����̘L�������f���郍�{�b�g���V�~�����[�g���A
�L���̒[�̃m�C�Y�I�h���g������l�ƃm�C�Y�̂���͈͑���l���g�p���܂��B���̗�ł́A�m�C�Y�̑���
�I�h���g���ƃZ���T�[�̓ǂݎ��l��Z�����邱�Ƃɂ��A�e�^�C���X�e�b�v�ł̃��{�b�g�̃|�[�Y��
�Ŗސ���iMLE�j���v�Z������@�������܂��B���{�b�g�͌��_����J�n���A�u-corridor_length�v�t���O��
�w�肳�ꂽ�Œ蒷�̃R���h�[�̏I�_�܂ňړ����܂��B��A�̃��[�V�����R�}���h�����s���āA
�u-pose_separation�v�t���O�Ŏw�肳�ꂽ�Œ蒷��O���Ɉړ����܂��B�ǂ̃|�[�Y�ŁA���ΓI��
���s��������l�ƁA�L���̒[�܂ł̋����͈̔͂̓ǂݎ��l���󂯎��܂��B�I�h���g���̓ǂݎ��l�́A
�u-odometry_stddev�v�t���O�Ŏw�肳�ꂽ�K�E�X�m�C�Y�ƕW���΍��ŕ`�悳��A�͈͂̓ǂݎ��l�́A
�u-range-stddev�v�t���O�Ŏw�肳�ꂽ�W���΍��œ��l�ɕ`�悳��܂��B���̖��ɂ�2�̃^�C�v��
�c��������܂��B
1�jOdometryConstraint�c���B����́A���{�b�g�̘A������|�[�Y����Ԃ̃I�h���g���ǂݎ��l��������܂��B
2�j�e�|�[�Y����ϑ����ꂽ�͈͂̓ǂݎ��l�̌덷���������RangeConstraint�c���B
OdometryConstraint�c�]�́A�Œ�p�����[�^�[�u���b�N�T�C�Y��1��
AutoDiffCostFunction�Ƃ��ă��f��������܂��B����́A��������鑊�΃I�h���g���ł��B
���{�b�g�̘A������|�[�Y�̃y�A�̊ԁB�ϑ����ꂽ���΃I�h���g���l�ƌv�Z���ꂽ���΃I�h���g���l�̍��́A
�I�h���g������l�̊��m�̕W���΍��ɂ���ďd�ݕt������ăy�i���e�B���ۂ����܂��B
RangeConstraint�c�]�́ADynamicAutoDiffCostFunction�Ƃ��ă��f��������܂��B
����́A���ΓI�ȑ��s�����̐���l�����v���āA���{�b�g�̐��肳�ꂽ�O���[�o���|�[�Y���v�Z���A
���ɁA�\�z�����͈͂̓ǂݎ��l���v�Z���܂��B���ɁA�ϑ����ꂽ�͈͂̓ǂݎ��l�Ɨ\�z���ꂽ
�͈͂̓ǂݎ��l�̍����A�Z���T�[�̓ǂݎ��l�̕W���΍��ɂ���ďd�ݕt������ăy�i���e�B��
�ۂ����܂��B���{�b�g�̃|�[�Y�̐��̓R���p�C�����ɕs���ł��邽�߁A���̃R�X�g�֐���
DynamicAutoDiffCostFunction�Ƃ��Ď�������܂��B��̏o�͂́A�I�h���g���Ɣ͈͂̓ǂݎ��l��
�����l�A����у��{�b�g�̂��ׂẴ|�[�Y�͈̔͂ƃI�h���g���̃G���[�ł��BMLE���v�Z������A
�v�Z���ꂽ�|�[�Y�ƏC�����ꂽ�I�h���g���l���A�Ή�����͈͂ƃI�h���g���G���[�ƂƂ��Ɉ������܂��B
�m�C�Y�̑����V�X�e����MLE�Ƃ��āA�G���[�̓[���Ɍ������܂��񂪁A�I�h���g���̐���l�́A
���{�b�g�̂��ׂẴI�h���g���Ɣ͈͂̓ǂݎ��l�̍Ŗޖ@���ő剻����悤�ɍX�V����邱�Ƃɒ��ӂ��Ă��������B

p_0�A..�Ap_N���iN + 1�j���{�b�g�̃|�[�Y�Ƃ��܂��B�����ŁA���{�b�g��p_0����n�܂�A
p_N�ŏI���L�����ړ����܂��Bp_0�����W�n�̌��_�ł���Ɖ��肵�܂��B�I�h���g��u_i�́A
�|�[�Yp_�ii-1�j��p_i�̊ԂŊϑ����ꂽ���΃I�h���g���ł���A�͈͓ǂݎ��ly_i�́A
�|�[�Yp_i����̃R���h�[�̒[�͈͓̔ǂݎ��l�ł��B�I�h���g���Ɣ͈͂̓ǂݎ��l��
�ǂ�����m�C�Y�������ł����A�M�O���œK�������悤�ɁA�C�����ꂽ�I�h���g���l
u * _0����u * _�iN-1�j�̍Ŗސ���iMLE�j���v�Z�������Ǝv���܂�
�FBelief�iu * _�i0�FN-1�j| u_�i0�FN-1�j�Ay_�i0�FN-1�j�j1�B
= P�iu * _�i0�FN-1�j| u_�i0�FN-1 �j�Ay_�i0�FN-1�j�j2�B
\ propto P�iy_�i0�FN-1�j| u * _�i0�FN-1�j�Au_�i0�FN-1�j�jP�iu * _�i0�FN-1�j| u_�i0�FN-1�j�j3�B
= \ prod_i {P�iy_i | u * _�i0�Fi�j�jP�iu * _i | u_i�j} 4.
�����ŁA���t�������u�i0�Fi�j�v�́A���̕ϐ��̂��ׂẴ^�C���X�e�b�v0����i�܂ł�
�G���g�����������߂̏ȗ��`�Ƃ��Ďg�p����܂��B�x�C�Y�̒藝�́A���𓱏o���邽�߂�
�g�p����܂��B2����3�ł���A�I�h���g���ϑ��Ɣ͈͓ǂݎ��̓Ɨ����́A3����4��
���o���邽�߂ɔr������܂��B���������āA�X�P�[���܂ł̐M�O�́A�e�|�[�Y��2�A
�e�|�[�Y��2�̗p��̐ςƂ��Ĉ�����������܂��B�p��͈͂̓ǂݎ��ɂ�1�̗p��
P�iy_i | u * _�i0�Fi�j�j�Ƒ��s��������̓ǂݎ��ɂ�1�̗p��P�iu * _i | u_i�j��
����܂��B�͈͂̓ǂݎ��̗p��͂Ɉˑ����邱�Ƃɒ��ӂ��Ă��������B���ׂĂ�
�I�h���g���lu * _�i0�Fi�j�A�I�h���g����P�iu * _i | u_i�j�́A�P��̒lu_i�݂̂Ɉˑ����܂��B
�͈͓ǂݎ��l�ƃI�h�G���g���[�m�����̗��������K���z�Ƃ��ă��f��������܂��B
�A����ь`���͎��̂Ƃ���ł��B
p�ix�j\ propto \ exp {-�i�ix --x_mean�j/ x_stddev�j^ 2}
�����ŁAx��MLE�I�h���g��u *�܂��͔͈͓ǂݎ��ly�̂����ꂩ���w���Ax_mean��
�I�h���g�����̑Ή����镽�ϒlu�ł��B �A�����y_expected�́A�ȑO�̂��ׂĂ�
�I�h���g�����Ɋ�Â��\�z�͈͂̓ǂݎ��l�ł��B���������āAMLE�́A�ŏ�������l
x *�������邱�Ƃɂ���Č������܂��B
x* = \ arg \ min {�i�ix --x_mean�j/ x_stddev�j^ 2}
�́A����`�ŏ����`���ł���ACeres�ɂ������ɓK���Ă��܂��B �B
����`�����́Ax_mean�̌v�Z���琶���܂��B
Ceres���œK������c���̎c���i�ix --x_mean�j/ x_stddev�j�B
�O�q�̂悤�ɁA�e�|�[�Y�̃I�h���g������1�̕ϐ��݂̂Ɉˑ����AAutoDiffCostFunction�ɂ����
�v�Z����܂����A�͈͓ǂݎ��̍��́A�ȑO�̂��ׂẴI�h���g���ϑ��Ɉˑ����܂��B
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

	// �Z���X���ɒǉ�����̂ɕ֗���RangeConstraint����CostFunction���쐬����t�@�N�g���[���\�b�h�D
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
		// �|�[�Yi����RangeConstraint��DynamicAutoDiffCostFunction���쐬���A�ǉ����܂��B
		vector<double*> param;
		Calculator::dyncostfunc* f =
			Calculator::Create(
				i, r_read[i], &o_val, &param);
		problem.AddResidualBlock(f, NULL, param);

		// �|�[�Y i �� OdometryConstraint �p�� AutoDiffCostFunction ���쐬���A�ǉ����܂��B
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

