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
using Eigen::IOFormat;//���l�m�F�p�D��Ő�Ώ���

class B4P_Calculator {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:
	static B4P_Bearing b4p;
	double dyn_h;		// �ŏ��X�e�b�v�T�C�Y
	double dyn_hmin;	// ���e�덷
	double dyn_ep;		// �������l
	double dyn_tr;		// �����X�e�b�v�T�C�Y
	int dyn_n;			// �ϐ�y�̗v�f��
	int ipar[128];		// intel OED solver �̐ݒ�l
	double *dpar;		// intel OED solver �̍�Ɣz��
	double t_step;		// ����̓T���v�����O���ԃX�e�b�v[s]
	double calctime;	// ����͌v�Z���[s]
	int nX;				// �ϐ� X �̔z�񒷂�

	struct Static {
		double eps[5];			// �É��臒l
		int iter1;				// �ő�J��Ԃ���
		int iter2;				// �ő厎�s��
		double rs;				// �����X�e�b�v�����D
		double jac_eps;			// ���R�r�A���v�Z�̐��x�D		
		int n;					// ���͔z��̒���
		int m;					// �o�͔z��̒���
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

