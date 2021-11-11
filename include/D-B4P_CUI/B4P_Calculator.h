#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Eigen\Dense>
#include<iostream>
#include<fstream>
#include "B4P_Bearing.h"
#include "B4P_FileIn.h"
#include "intel_ode.h"
#include "mkl_rci.h"

using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::IOFormat;//数値確認用．後で絶対消す

class B4P_Calculator {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:
	static B4P_Bearing b4p;
	double dyn_h;		// 最小ステップサイズ
	double dyn_hmin;	// 許容誤差
	double dyn_ep;		// しきい値
	double dyn_tr;		// 初期ステップサイズ
	int dyn_n;			// 変数yの要素長
	int ipar[128];		// intel OED solver の設定値
	double *dpar;		// intel OED solver の作業配列
	double t_step;		// 動解析サンプリング時間ステップ[s]
	double calctime;	// 動解析計算区間[s]
	int nX;				// 変数 X の配列長さ

	struct Static {
		double eps[5];			// 静解析閾値
		int iter1;				// 最大繰り返し数
		int iter2;				// 最大試行回数
		double rs;				// 初期ステップ条件．
		double jac_eps;			// ヤコビアン計算の精度．		
		int n;					// 入力配列の長さ
		int m;					// 出力配列の長さ
	} SttSet[3];

private:
	static void Dyn_Eq(int *n, double *t0, double *y, double *dydt);
	static void Dyn_void(void);
	
public:
	int Stt_solve_stf(double*x);
	void Stt_init(const B4P_FileIn::Static * SttSet);
	static void Stt_Eq_stf(int *m, int *n, double *x, double *f);
	void Dyn_init(const B4P_FileIn::Dynamic & DynSet);
	void Dyn_solve(double*x0, double*t);
	~B4P_Calculator();
};

