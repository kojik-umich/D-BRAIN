#pragma once
#include <direct.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdlib.h>

#include "B4P_In.h"
#include "Numeric.h"
#include "FileIn.h"
#include "Unit.h"

using std::ofstream;

class B4P_FileIn : public B4P_In {
public:
	struct FileName {
		string Temp;
		string Ball;
		string Inner;
		string Outer;
		string Cage;
		string BallInnerPair;
		string BallOuterPair;
		string BallCagePair;
	} FN;

	struct Dynamic {
		int sampling;			// 動解析結果出力サンプリング数
		double calctime;		// 動解析計算時間（計算開始時点の時刻を0sとおいた時の終了時刻）[s]
		double h;				// 動解析初期時間ステップサイズ[s]
		double hmin;			// 動解析最小時間ステップサイズ[s]
		double ep;				// 動解析相対許容誤差
		double tr;				// 動解析相対許容誤差しきい値
		double x[3];			// 内輪変位 [mm]
		double ax[3];			// 内輪軸方向 [-]
		bool v_is_Locked[3];	// 並進方向拘束条件
		bool w_is_Locked[3];	// 回転方向拘束条件
		int nX;					// 変数 x の長さ
		bool fromcontinuation;	// 続き計算
		bool stopcalculation;	// 計算打ち切り
		double dTierr;			// 許容トルク誤差（計算打ち切り設定時）
		int stp;				// 打ち切りステップ数
	} DynSet;

	struct Static {
		double eps[5];			// 静解析閾値
		int iter1;				// 最大繰り返し数
		int iter2;				// 最大試行回数
		double rs;				// 初期ステップ条件．
		double jac_eps;			// ヤコビアン計算の精度．		
		int nX;					// 変数 x の長さ
	} SttSet[3];

	struct OutputSlice {
		int n;					// スライス荷重を出力する玉数
		int *list;				// スライス荷重を出力する玉番号
	} OutSlice;
	


public:
	void Nothing(void);		// B4P_Inの具象クラスであることを示すため，意味のない関数を定義する．

	bool read_input_all(char fpath_d4bin_csv[]);

	bool read_BallNum(const vector<vector<string>> &inp_data);
	bool read_SetCage(const vector<vector<string>> &inp_data);
	bool read_Ball(const vector<vector<string>> &inp_data, int i);
	bool read_LoadIn(const vector<vector<string>> &inp_data);

	bool read_Inner(const vector<vector<string>> &inp_data);
	bool read_Outer(const vector<vector<string>> &inp_data);

	bool read_Oil(const vector<vector<string>> &inp_data);

	bool read_SnapCage(const vector<vector<string>> &inp_data);
	bool read_ContaBI(const vector<vector<string>> &inp_data);
	bool read_ContaBO(const vector<vector<string>> &inp_data);
	bool read_ContaBC(const vector<vector<string>> &inp_data);
	bool read_ContaCI(const vector<vector<string>> &inp_data);
	bool read_ContaCO(const vector<vector<string>> &inp_data);

	bool read_Gravity(const vector<vector<string>> &inp_data);
	bool read_Rotation(const vector<vector<string>> &inp_data);
	bool read_DynSet(const vector<vector<string>> &inp_data);
	bool read_SttSet(const vector<vector<string>> &inp_data, int i);

	bool read_InnnerBound(const vector<vector<string>> &inp_data);
	bool read_InnnerPosition(const vector<vector<string>> &inp_data);

	bool read_FileName(const vector<vector<string>> &inp_data, string fname_out);
	bool read_RollingResistance(const vector<vector<string>> &inp_data);
	bool read_Coulomb(const vector<vector<string>>&inp_data);
	bool read_Hysteresis(const vector<vector<string>>&inp_data);
	bool read_FilmThickness(const vector<vector<string>>&inp_data);
	bool read_Dimension(const vector<vector<string>> &inp_data);
	bool read_FromContinuation(const vector<vector<string>> &inp_data);
	bool read_StopCalculation(const vector<vector<string>> &inp_data);
	bool read_OutputSlice(const vector<vector<string>> &inp_data);
	static double calc_cos_alp0(double R, double D, double RO_x);
};

//enum AnalysisMode {
//	Output_only = 0,			// 出力のみ
//	Dyn_from_beginning = 1,		// 動解析（初めから計算）
//	Dyn_from_continuation = 2,	// 動解析（続きから計算）
//}AnaMode;

