#pragma once
#include "pch.h"

class B4P_StabBallInnerRingPair : public B4P_BallInnerRingPair {
public:
	virtual void calc_force(Vector3d&Fbi, Vector3d&Nbi, Vector3d&Fib, Vector3d&Nib) {
		Vector3d ones = Vector3d::Ones();
		Fbi = 1.0 * ones;
		Nbi = 10.0 * ones;
		Fib = 100.0 * ones;
		Nib = 1000.0 * ones;
		return;
	};
};

class B4P_StabBallOuterRingPair : public B4P_BallOuterRingPair {
public:
	virtual void calc_force(Vector3d&Fbi, Vector3d&Nbi, Vector3d&Fib, Vector3d&Nib) {
		Vector3d ones = Vector3d::Ones();
		Fbi = 1.0 * ones;
		Nbi = 10.0 * ones;
		Fib = 100.0 * ones;
		Nib = 1000.0 * ones;
		return;
	};
};

class B4P_StabBallCagePair : public B4P_BallCagePair {
public:
	virtual void calc_force(Vector3d&Fbi, Vector3d&Nbi, Vector3d&Fib, Vector3d&Nib) {
		Vector3d ones = Vector3d::Ones();
		Fbi = 1.0 * ones;
		Nbi = 10.0 * ones;
		Fib = 100.0 * ones;
		Nib = 1000.0 * ones;
		return;
	};
};







