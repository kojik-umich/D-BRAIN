#pragma once

#include <iostream>
#include <fstream>

#include "Unit.h"
#include "FileOut.h"

#include "BS_In.h"
#include "BS_Out.h"
#include "BS_FileIn.h"

using std::string;
using std::to_string;

class BS_FileOut : public BS_Out {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:
	void Nothing(void);		// BS_Inの具象クラスであることを示すため，意味のない関数を定義する．

private:
	int np;

	string _01_Ball_csv;
	string _02_Nut_csv;
	string _03_Shaft_csv;
	string _04_BallNutPair_csv;
	string _05_BallShaftPair_csv;
	string _06_BallBallPair_csv;

	VectorXd BL_;
	VectorXd NT_;
	VectorXd ST_;
	int ncn;
	VectorXd NT_CY_;
	int ncs;
	VectorXd ST_CY_;
	VectorXd BNP_;
	VectorXd BSP_;
	VectorXd BBP_;

	VectorXd BL_All;
	VectorXd NT_All;
	VectorXd ST_All;
	VectorXd BNP_All;
	VectorXd BSP_All;
	VectorXd BBP_All;

public:
	void init(const BS_FileIn&FI);

	void write_AllHeader(void);
	void write_AllParams(double Time);
	static void read_Params(string name, int nX, double &t, double *Params);

public:
	static int write_precision;
	static int iRG;
	static int iCY;
	static int iBCP;
	static int iBBP;

private:
	static VectorXd form_RG(const BS_Out::Rigid & RG);
	static VectorXd form_CY(const BS_Out::Cylinder & CY);
	static VectorXd form_BCP(const BS_Out::BallCylinderPair&BCP);
	static VectorXd form_BBP(const BS_Out::BallBallPair&BBP);
};

