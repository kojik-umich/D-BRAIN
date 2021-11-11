#pragma once

#include <iostream>
#include <fstream>

#include "Unit.h"
#include "FileOut.h"

#include "B4P_Out.h"
#include "B4P_FileIn.h"
#include "bal_Cage.h"	// 接触点数取得のため．

using std::string;
using std::to_string;

class B4P_FileOut : public B4P_Out {
public:
	void Nothing(void);		// B4P_Inの具象クラスであることを示すため，意味のない関数を定義する．

private:
	int msmax;
	int write_precision;
	int OutSlice_n;
	int *OutSlices;
	string _01_Ball_csv;
	string _02_Outer_csv;
	string _03_Inner_csv;
	string _04_Cage_csv;
	string _05_BallOuterPair_csv;
	string _06_BallInnerPair_csv;
	string _07_BallCagePair_csv;

	VectorXd Rigid_;
	VectorXd BRP_;
	VectorXd BCP_;

	VectorXd Balls_;
	VectorXd BRPs_;
	VectorXd BCPs_;
	VectorXd Ring_;
	VectorXd Slice_;

	static int Z;
	static int Rigid_size;
	static int BRP_size;
	static int BCP_size;
	static int Slice_size;

public:
	void init(const B4P_FileIn&FI);

	void write_AllHeader(void);
	void write_AllParams(double Time);
	string make_slices_header(string bi, string gj);

	void write_Header(string name);
	void write_Params(string name, double Time, const double*Params);
	void read_Params(string name, double& t, double*Params);
	static VectorXd form_RG(const B4P_Out::Rigid & RG);
	static VectorXd form_BRP(const B4P_Out::BallRingPair & BRP);
	static VectorXd form_BCP(const B4P_Out::BallCagePair & BCP);
	VectorXd form_SL(const B4P_Out::BallRingPair & BRP);
};

