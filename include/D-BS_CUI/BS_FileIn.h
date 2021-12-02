#pragma once
#include <sstream>		// String
#include <vector>		// vector
#include <map>			// Map

#include "BS_In.h"		// 継承元
#include "Numeric.h"	// 2乗の計算など
#include "FileIn.h"		// ファイル読み込み
#include "Unit.h"		// 単位変換
#include "Tribology.h"	// 換算質量など

using namespace std;

class BS_FileIn : public BS_In {
public:

	enum Mass {
		Reduced = 0,	// 換算質量で計算．
		Nut = 1,		// ナットの質量で計算
		Shaft = 2		// シャフトの質量で計算
	} mass;

	struct Initial {
		enum Preset {			// 初期値設定
			ReadPos		= 0,	// $$Positionから読み取り
			ReadTemp	= 1		// Tempファイルから読み取り
		}preset;
		double x0[3];	// [m/s]:	初期入力変位
		double ax0[3];	// [rad/s]:	初期入力姿勢
	} initial;

	struct Preload {
		enum Mode {
			distance = 0,
			load = 1
		} mode;
	}preload;

	struct Output {
		string _01;
		string _02;
		string _03;
		string _04;
		string _05;
		string _06;
		string temp;
		bool deletelastout; // 前回計算結果を消去
	}output;

	bool runStatic;	// 静解析の実行の有無
	struct Static {
		bool run[5];		// 静解析の実行の有無
		struct Set {
			int n;				// 入力配列サイズ．
			int m;				// 出力配列サイズ．
			int iter1;			// 最大繰り返し数．
			int iter2;			// 最大試行回数．
			double rs;			// 初期ステップ条件．
			double jac_eps;		// ヤコビアン計算の精度．
			double eps[6];		// ソルバの設定（長さ6の配列）．各要素の意味はマニュアル参照．
		} set[3];
		double v0;				// シャフト進行速度[m/s]
		double w0;				// シャフト回転速度[rad/s]
		double wn;				// ナット回転速度[rad/s]
	} stt;

	struct Dynamic {
		int wxt_n, vxd_n, ft_n;			// 入力パラメータの個数(wxt, vxtの読み込みにしか使わない)
		map<double, double>wxt;			// 角速度時間変化(first: 時間[s], second:角速度[rad/s])
		map<double, double>vxt;			// 速度時間変化　(first: 時間[s], second:速度[m/s])
		map<double, vector<double>>ft;	// 荷重時間変化　(first: 時間[s], second:荷重[N]/[Nm])
		
		struct Set {
			double calctime;	// 動解析計算時間（計算開始時点の時刻を0sとおいた時の終了時刻）[s]
			int    sampling;	// 動解析結果出力サンプリング数
			double h;			// 動解析初期時間ステップサイズ[s]
			double hmin;		// 動解析最小時間ステップサイズ[s]
			double ep;			// 動解析相対許容誤差
			double tr;			// 動解析相対許容誤差しきい値

			bool stopcalc;		// 計算自動停止する(1)，しない(0)
			double dTerr;		// 自動停止する場合のトルク誤差閾値[Nm]
			int stp;			// 何step連続で閾値を下回ったら計算終了するか
		} set[2];
		
	} dyn;

public:
	void Nothing(void);		// BS_Inの具象クラスであることを示すため，意味のない関数を定義する．

public:
	void read_input_all(const char fpath_d4bin_csv[]);

private:
	bool read_allocate1(const vector<vector<string>>& inp_data);
	bool read_allocate2(const vector<vector<string>>& inp_data);
	bool read_allocate3(const vector<vector<string>>& inp_data);
	bool read_NutNum(const vector<vector<string>>& inp_data, int & nutnum);
	bool read_SpiralNum(const vector<vector<string>>& inp_data, int & spiralnum);
	bool read_CircuitNum(const vector<vector<string>>& inp_data, int & circuitnum);
	bool read_Circuit(const vector<vector<string>>& inp_data, int i, int & ballnum);
	bool read_PreLoad(const vector<vector<string>>& inp_data);
	bool read_RollingResistance(const vector<vector<string>>&inp_data);
	bool read_Coulomb(const vector<vector<string>>& inp_data);
	bool read_Ellipse(const vector<vector<string>>& inp_data);
	bool read_FilmThickness(const vector<vector<string>>& inp_data);
	bool read_Hysteresis(const vector<vector<string>>& inp_data);
	bool read_Dimension(const vector<vector<string>>&inp_data);
	bool read_Ball(const vector<vector<string>>& inp_data);
	bool read_Shaft(const vector<vector<string>>& inp_data, double & Shaft_PCD);
	bool read_ShaftSpiral(const vector<vector<string>>& inp_data, int i, double PCD);
	bool read_BallShaftPair(const vector<vector<string>>& inp_data);
	bool read_Nut(const vector<vector<string>>& inp_data, int i, double & Nut_PCD);
	static double calc_Cylinder_m(double ri, double ro, double l, double rho);
	static double calc_Cylinder_Ix(double ri, double ro, double l, double rho);
	static double calc_Cylinder_Iyz(double ri, double ro, double l, double rho);
	static double calc_Cylinder_Iyz_l(double m, double x, double l);
	bool read_NutSpiral(const vector<vector<string>>& inp_data, int i, int j, double PCD);
	bool read_BallNutPair(const vector<vector<string>>& inp_data);
	bool read_ShaftMassSet(const vector<vector<string>>& inp_data);
	bool read_Oil(const vector<vector<string>>& inp_data);
	bool read_Gravity(const vector<vector<string>>& inp_data);
	bool read_PositionSet(const vector<vector<string>>& inp_data);
	bool read_Position(const vector<vector<string>>& inp_data);
	bool read_SttLoadNum(const vector<vector<string>>& inp_data, int & sttloadnum);
	bool read_SttLoad(const vector<vector<string>>& inp_data, int i);
	bool read_SttRotation(const vector<vector<string>>& inp_data);
	bool read_SttMode(const vector<vector<string>>& inp_data);
	bool read_SttSet(const vector<vector<string>>& inp_data, int i);
	bool read_DynSet(const vector<vector<string>>& inp_data, int i);
	bool read_Bound(const vector<vector<string>>& inp_data);
	void calc_Mass(void);
	bool read_Output(const vector<vector<string>>& inp_data);
	bool read_DynRotationStep(const vector<vector<string>>& inp_data);
	bool read_DynRotation(const vector<vector<string>>& inp_data, int i);
	bool read_DynLoadStep(const vector<vector<string>>& inp_data);
	bool read_DynLoad(const vector<vector<string>>& inp_data, int i);
	bool read_StopDynCalc(const vector<vector<string>>& inp_data, int i);
	static double calc_phase(double l, double x);
	void set_Circuitphase(int i);
	bool read_DeleteLastOutput(const vector<vector<string>>&inp_data);
	bool read_DynVelocityStep(const vector<vector<string>>&inp_data);
	bool read_DynVelocity(const vector<vector<string>>&inp_data, int i);
};
