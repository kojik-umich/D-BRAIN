#pragma once
#include "Numeric.h"
#include <Eigen\Dense>
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Quaterniond;

/// <summary>����</summary>
/// <remarks>�ψʂ�p���C���ʂ������C�x�N�g���̕ϊ���^���������ɂ�鎞�Ԕ����l���v�Z����D</remarks>
class Rigid {
public:
	Vector3d x;		// �ψ� [m]�i0:x,1:y,2:z�j
	Vector3d v;		// ���x [m/s]
	Quaterniond q;	// �p����\���N�H�[�^�j�I�� [-]
	Vector3d w;		// �p���x [rad/s]

	double   m;		// ����[kg]
	double   m_inv;	// ���ʂ̋t���D
	Vector3d I;		// �������[�����g[kg*m^2]
	Vector3d I_inv;	// �������[�����g�̋t��
	double   rho;	// ���x [kg/m^3]
	double   E;		// �����O�� [Pa]
	double   nu;	// �|�A�\���� [-]

	Vector3d ax0;	// �����̎������x�N�g���D(=[1,0,0])
	Vector3d ay0;	// ���ʊp180���̌a�����x�N�g���D����D(=[0,1,0])
	Vector3d az0;	// ���ʊp270���̌a�����x�N�g���D����D(=[0,0,1])

	Vector3d F;		// �󂯂Ă���O���׏d[N]�D�������W�n�D
	Vector3d T;		// �󂯂Ă���O���g���N[Nm]�D�������W�n�D

	//bool v_const;	// ���x�Œ�
	//bool w_const;	// �p���x�Œ�


	Eigen::Array<bool, 3, 1>x_const;		//	�S�������ix, y, z�����j
	Eigen::Array<bool, 3, 1>Rx_const;		//	�S�������iRx, Ry, Rz�����j

	static double l;	// �P�ʒ���[m]
	static double t;	// �P�ʎ���[s]
	static Vector3d g;	// �d�͉����x[m/ss]

public:
	Rigid(void);
	~Rigid(void);
	virtual void Nothing(void)=0;	// ���ۃN���X�ł��邱�Ƃ��������߁C�Ӗ��̂Ȃ��֐����`����D
	void set_param(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w);
	void get_param(Vector3d&x, Vector3d&v, Quaterniond&q, Vector3d&w);
	void set_y(const Vector3d & x, const Vector3d & v, const Quaterniond & q, const Vector3d & w);
	void get_y(Vector3d & x, Vector3d & v, Quaterniond & q, Vector3d & w);
	Vector3d get_dxdt(void);
	Vector3d get_dvdt(const Vector3d&F);
	Quaterniond get_dqdt(void);
	Vector3d get_dwdt(const Vector3d&T);
	void get_dydt(const Vector3d&F, const Vector3d&T, double*dydt);
	void set_const(const bool (&v_const)[3] ,const bool (&w_const)[3]);
	Vector3d to_mycoord(const Vector3d&x);
	Vector3d to_inecoord(const Vector3d&x);
	Vector3d to_myvelocity(const Vector3d&v);
	Vector3d to_inevelocity(const Vector3d&V);
	Vector3d to_myvector(const Vector3d&x);
	Vector3d to_inevector(const Vector3d&x);
	Vector3d surface_velocity(const Vector3d&x);
	Vector3d calc_Torque(const Vector3d&x, const Vector3d&F);
	Vector3d calc_TorqueDirection(const Vector3d & x, const Vector3d & u);
	Vector3d get_ax(void);
	Vector3d get_ay(void);
	Vector3d get_az(void);
	Vector3d get_mg(void);
	Vector3d remove_ax(const Vector3d&x);
	Vector3d remove_vector(const Vector3d&x);
	void set_ax(const Vector3d&ax);
	void set_FT(const Vector3d&F, const Vector3d&T);
	void set_mI(double m, const Vector3d&I);
	void save(double * x, double * v, double * q, double * w, double * ax, double * F, double * T);
};

