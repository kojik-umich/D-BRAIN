#pragma once
#include <iostream>
#include <Eigen\Dense>
#include "Ball.h"
#include "Tribology.h"
#include "bal_Cage.h"
#include "B4P_In.h"
#include "B4P_Out.h"
using Eigen::Vector3d;
using Eigen::ColPivHouseholderQR;

// �{�[��-�ێ���Ԃ̃��[�v�s�ϗʂ�ۑ����Ă����N���X�DCG=cage�D
class B4P_BallCagePair{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
private:
	Tribology::DampingForce*DF;

	Ball       *BL;
	bal_Cage   *CG;
	int         np;		// �{�[�����Ή�����ێ���|�P�b�g�ԍ��D
	double      mu;		// �ێ���-�{�[���̃N�[�������C�W���D(typical:0.1~0.3)
	double      zeta;	// �ێ���-�{�[���̃_���s���O�W���D(typical:~0.2)
	double      m;		// ���Z���ʁD
	Vector3d    p[_MAX_CONTACT_];	// �ڐG�_�ʒu�i�������W�n�j
	Vector3d	Fn[_MAX_CONTACT_];	// �ڐG�_�׏d�i�������W�n�j
	Vector3d	Fs[_MAX_CONTACT_];	// �ڐG�_���C�i�������W�n�j
	double		dx[_MAX_CONTACT_];	// �ڋߗ�[m]
	int			ptt[_MAX_CONTACT_];	// �ڐG�p�^�[���i�{����enum�ŊǗ����ׂ������C�ێ���ɂ���Ē�`���قȂ邽�߁Cint�ŊǗ��j


public:
	void init(const B4P_In&FI);
	void link(Ball*BL, bal_Cage*CG, int np);
	virtual void calc_force(Vector3d&Fbc,Vector3d&Tbc,Vector3d&Fcb,Vector3d&Tcb);
	void save(B4P_Out::BallCagePair&BCP);
	Vector3d get_us(const Vector3d&p);
};

