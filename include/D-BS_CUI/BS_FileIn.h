#pragma once
#include <sstream>		// String
#include <vector>		// vector
#include <map>			// Map

#include "BS_In.h"		// �p����
#include "Numeric.h"	// 2��̌v�Z�Ȃ�
#include "FileIn.h"		// �t�@�C���ǂݍ���
#include "Unit.h"		// �P�ʕϊ�
#include "Tribology.h"	// ���Z���ʂȂ�

using namespace std;

class BS_FileIn : public BS_In {
public:

	enum Mass {
		Reduced = 0,	// ���Z���ʂŌv�Z�D
		Nut = 1,		// �i�b�g�̎��ʂŌv�Z
		Shaft = 2		// �V���t�g�̎��ʂŌv�Z
	} mass;

	struct Initial {
		enum Preset {			// �����l�ݒ�
			ReadPos		= 0,	// $$Position����ǂݎ��
			ReadTemp	= 1		// Temp�t�@�C������ǂݎ��
		}preset;
		double x0[3];	// [m/s]:	�������͕ψ�
		double ax0[3];	// [rad/s]:	�������͎p��
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
		bool deletelastout; // �O��v�Z���ʂ�����
	}output;

	bool runStatic;	// �É�͂̎��s�̗L��
	struct Static {
		bool run[5];		// �É�͂̎��s�̗L��
		struct Set {
			int n;				// ���͔z��T�C�Y�D
			int m;				// �o�͔z��T�C�Y�D
			int iter1;			// �ő�J��Ԃ����D
			int iter2;			// �ő厎�s�񐔁D
			double rs;			// �����X�e�b�v�����D
			double jac_eps;		// ���R�r�A���v�Z�̐��x�D
			double eps[6];		// �\���o�̐ݒ�i����6�̔z��j�D�e�v�f�̈Ӗ��̓}�j���A���Q�ƁD
		} set[3];
		double v0;				// �V���t�g�i�s���x[m/s]
		double w0;				// �V���t�g��]���x[rad/s]
		double wn;				// �i�b�g��]���x[rad/s]
	} stt;

	struct Dynamic {
		int wxt_n, vxd_n, ft_n;			// ���̓p�����[�^�̌�(wxt, vxt�̓ǂݍ��݂ɂ����g��Ȃ�)
		map<double, double>wxt;			// �p���x���ԕω�(first: ����[s], second:�p���x[rad/s])
		map<double, double>vxt;			// ���x���ԕω��@(first: ����[s], second:���x[m/s])
		map<double, vector<double>>ft;	// �׏d���ԕω��@(first: ����[s], second:�׏d[N]/[Nm])
		
		struct Set {
			double calctime;	// ����͌v�Z���ԁi�v�Z�J�n���_�̎�����0s�Ƃ��������̏I�������j[s]
			int    sampling;	// ����͌��ʏo�̓T���v�����O��
			double h;			// ����͏������ԃX�e�b�v�T�C�Y[s]
			double hmin;		// ����͍ŏ����ԃX�e�b�v�T�C�Y[s]
			double ep;			// ����͑��΋��e�덷
			double tr;			// ����͑��΋��e�덷�������l

			bool stopcalc;		// �v�Z������~����(1)�C���Ȃ�(0)
			double dTerr;		// ������~����ꍇ�̃g���N�덷臒l[Nm]
			int stp;			// ��step�A����臒l�����������v�Z�I�����邩
		} set[2];
		
	} dyn;

public:
	void Nothing(void);		// BS_In�̋�ۃN���X�ł��邱�Ƃ��������߁C�Ӗ��̂Ȃ��֐����`����D

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
