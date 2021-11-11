#pragma once
#include "Rigid.h"
#include "Ball.h"
#include "BS_In.h"
#include "BS_Cylinder.h"
#include <vector>

class BS_Nut : public Rigid {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:
	int nCY;			// num-cylinder．シングルで１，ダブルで２．
	BS_Cylinder*CY;		// 螺旋溝円筒部材．

public:
	virtual void Nothing(void) = 0;	// まだ抽象クラスであるため，実態は定義しない．
	virtual void init(const std::vector<BS_In::Cylinder>& nut, double w0) = 0;
	void allocate(const std::vector<BS_In::Cylinder>& cylinder);
	virtual void set_y_(const Vector3d & x, const Vector3d & v, const Quaterniond & q, const Vector3d & w) = 0;

	BS_Nut();
	~BS_Nut();
};

class BS_SingleNut : public BS_Nut {

public:
	void Nothing(void);	// 具象クラスであることを示すため，意味のない関数を定義する．
	void init(const std::vector<BS_In::Cylinder>& nut, double w0);
	void set_y_(const Vector3d & x, const Vector3d & v, const Quaterniond & q, const Vector3d & w);
};


class BS_DoubleNut : public BS_Nut {
public:
	void Nothing(void);	// 具象クラスであることを示すため，意味のない関数を定義する．
	void init(const std::vector<BS_In::Cylinder>& nut, double w0);

private:
	Vector3d dx[2];		// 与圧距離[m]．シングルナットでは使用しない．
	void set_dx(void);
	void set_y_(const Vector3d & x, const Vector3d & v, const Quaterniond & q, const Vector3d & w);
};

