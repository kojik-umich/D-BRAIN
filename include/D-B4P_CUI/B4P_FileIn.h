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
		int sampling;			// ����͌��ʏo�̓T���v�����O��
		double calctime;		// ����͌v�Z���ԁi�v�Z�J�n���_�̎�����0s�Ƃ��������̏I�������j[s]
		double h;				// ����͏������ԃX�e�b�v�T�C�Y[s]
		double hmin;			// ����͍ŏ����ԃX�e�b�v�T�C�Y[s]
		double ep;				// ����͑��΋��e�덷
		double tr;				// ����͑��΋��e�덷�������l
		double x[3];			// ���֕ψ� [mm]
		double ax[3];			// ���֎����� [-]
		bool v_is_Locked[3];	// ���i�����S������
		bool w_is_Locked[3];	// ��]�����S������
		int nX;					// �ϐ� x �̒���
		bool fromcontinuation;	// �����v�Z
		bool stopcalculation;	// �v�Z�ł��؂�
		double dTierr;			// ���e�g���N�덷�i�v�Z�ł��؂�ݒ莞�j
		int stp;				// �ł��؂�X�e�b�v��
	} DynSet;

	struct Static {
		double eps[5];			// �É��臒l
		int iter1;				// �ő�J��Ԃ���
		int iter2;				// �ő厎�s��
		double rs;				// �����X�e�b�v�����D
		double jac_eps;			// ���R�r�A���v�Z�̐��x�D		
		int nX;					// �ϐ� x �̒���
	} SttSet[3];

	struct OutputSlice {
		int n;					// �X���C�X�׏d���o�͂���ʐ�
		int *list;				// �X���C�X�׏d���o�͂���ʔԍ�
	} OutSlice;
	


public:
	void Nothing(void);		// B4P_In�̋�ۃN���X�ł��邱�Ƃ��������߁C�Ӗ��̂Ȃ��֐����`����D

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
//	Output_only = 0,			// �o�͂̂�
//	Dyn_from_beginning = 1,		// ����́i���߂���v�Z�j
//	Dyn_from_continuation = 2,	// ����́i��������v�Z�j
//}AnaMode;

