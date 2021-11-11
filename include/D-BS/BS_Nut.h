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
	int nCY;			// num-cylinder�D�V���O���łP�C�_�u���łQ�D
	BS_Cylinder*CY;		// �����a�~�����ށD

public:
	virtual void Nothing(void) = 0;	// �܂����ۃN���X�ł��邽�߁C���Ԃ͒�`���Ȃ��D
	virtual void init(const std::vector<BS_In::Cylinder>& nut, double w0) = 0;
	void allocate(const std::vector<BS_In::Cylinder>& cylinder);
	virtual void set_y_(const Vector3d & x, const Vector3d & v, const Quaterniond & q, const Vector3d & w) = 0;

	BS_Nut();
	~BS_Nut();
};

class BS_SingleNut : public BS_Nut {

public:
	void Nothing(void);	// ��ۃN���X�ł��邱�Ƃ��������߁C�Ӗ��̂Ȃ��֐����`����D
	void init(const std::vector<BS_In::Cylinder>& nut, double w0);
	void set_y_(const Vector3d & x, const Vector3d & v, const Quaterniond & q, const Vector3d & w);
};


class BS_DoubleNut : public BS_Nut {
public:
	void Nothing(void);	// ��ۃN���X�ł��邱�Ƃ��������߁C�Ӗ��̂Ȃ��֐����`����D
	void init(const std::vector<BS_In::Cylinder>& nut, double w0);

private:
	Vector3d dx[2];		// �^������[m]�D�V���O���i�b�g�ł͎g�p���Ȃ��D
	void set_dx(void);
	void set_y_(const Vector3d & x, const Vector3d & v, const Quaterniond & q, const Vector3d & w);
};

