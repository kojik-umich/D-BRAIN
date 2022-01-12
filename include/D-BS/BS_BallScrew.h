#pragma once
#include <Eigen\Dense>
#include "BS_In.h"
#include "BS_Out.h"
#include "Ball.h"
#include "BS_Circuit.h"
#include "BS_Nut.h"
#include "BS_Shaft.h"
#include "BS_BallCylinderPair.h"

class BS_BallScrew {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
	int nCC;				// 回路の数

	BS_Circuit *CC;			// 回路
	BS_Nut  *NT;			// ナット（シングル/ダブル対応のため，ポインタで所持しています．）
	BS_Shaft ST;			// シャフト（ねじ軸）

	int nP;					// num-pair．全ての回路に含まれる玉数
	BS_BallNutPair*BNP;		// 玉-ナット間接触
	BS_BallShaftPair*BSP;	// 玉-シャフト間接触
	BallBallPair*BBP;		// 玉-玉間接触（未実装）

public:
	struct Load {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
		Vector3d F;			// 外部荷重 [N]
		Vector3d T;			// 外部トルク [Nm]
	} LD;

public:
	void allocate(const BS_In & IN);
	void link(const BS_In & IN);
	void init_position(double wn, double ws);
	void init_Load(const BS_In & IN);
	void init(const BS_In & IN, double v0, double w0, double wn);
	void preset_y0(double dx0, double dth0, double dx1);
	void pure_Rolling(void);
	bool get_y1(int ib, double * y2);
	void set_y1(int ib, const double * y1);
	void get_F1(int ib, double * f1);
	void lock_y0(const double * x0, const double * ax0, double v0, double w0);
	void preset_y0_F(double dx0);
	void preset_y0_T(double dth0);
	void preset_y0_x(double dx);
	void get_y0(double*y0);
	void set_y0(const double * y0, double v0, double w0);
	void get_F0(double*f0);

	void get_F1(double*f1);

	void get_y2(double*y0);
	void set_y2(const double * y2, double v0, double w0);
	void get_F2(double*f2);

	void set_dyn_y0(const double * y0);
	void init_dyn0(void);
	void deinit_dyn0(double v0, double w0);
	void get_dyn_y0(double * y0);
	void get_dyn_dydt0(double * y);

	void set_dyn_y1(const double * y);
	void get_dyn_y1(double * y);
	void get_dyn_dydt1(double * dydt, double v, double w);
	//void Shaft_Lock(void);
	void save(BS_Out&OUT);
	void set_load(double *F, double *T);
	BS_BallScrew();
	~BS_BallScrew();



private:
	static double Hertz(double k, double dx) {
		return k * (
			(dx > 0)
			? std::pow(dx, 1.5)
			: std::pow(-dx, 1.5) * -1e-4
			);
	};
	const double k0 = 1.0;
	const double k1 = 2.0;
	const double k2 = 4.0;
	const double F0 = 1.0;

	double x0;
	double x1;

public:

	void getPosition(std::vector<double>& x) const {
		x[0] = this->x0;
		x[1] = this->x1;
	};

	void setPosition(double const* const x) {
		this->x0 = x[0];
		this->x1 = x[1];
	};

	void getForce(double * F) {
		double dx0 = this->x0;
		double dx1 = this->x0 - this->x1;
		double dx2 = this->x1;

		F[0] = -Hertz(this->k0, dx0) - Hertz(this->k1, dx1) + this->F0;
		F[1] = Hertz(this->k1, dx1) - Hertz(this->k2, dx2);
	};
};






