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

class BS_SimulinkIn : public BS_In {

public:
	void Nothing(void);		// BS_Inの具象クラスであることを示すため，意味のない関数を定義する．

public:
	void read_csv(const char fpath_d4bin_csv[]);

private:
	bool read_allocate1(const vector<vector<string>>& inp_data);
	bool read_allocate2(const vector<vector<string>>& inp_data);
	bool read_NutNum(const vector<vector<string>>& inp_data, int & nutnum);
	bool read_SpiralNum(const vector<vector<string>>& inp_data, int & spiralnum);
	bool read_CircuitNum(const vector<vector<string>>& inp_data, int & circuitnum);
	bool read_Circuit(const vector<vector<string>>& inp_data, int i, int & ballnum);
	bool read_RollingResistance(const vector<vector<string>>&inp_data);
	bool read_Coulomb(const vector<vector<string>>& inp_data);
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
	bool read_Oil(const vector<vector<string>>& inp_data);
	bool read_Gravity(const vector<vector<string>>& inp_data);
	static double calc_phase(double l, double x);
	bool read_PreLoad(const vector<vector<string>>& inp_data);
	void set_Circuitphase(int i);
};
