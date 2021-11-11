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

	int i2;				// Step2で計算対象とする玉番号






	// 外部荷重
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
	void preset_y0(double dx0, double dth0);
	void lock_y0(const double * x0, const double * ax0, double v0, double w0);
	void preset_y0_F(double dx0);
	void preset_y0_T(double dth0);
	void preset_y0_x(double dx);
	void get_y0(double*y0);
	void set_y0(const double * y0, double v0, double w0);
	void get_F0(double*f0);

	void get_F1(double*f1);

	void set_i2(int i2);
	void get_y2(double*y0);
	void set_y2(const double*y0);
	void get_F2(double*f2);


	void set_y(const double * y);
	void get_y(double * y);
	void get_dydt(double * dydt, double v, double w);
	//void Shaft_Lock(void);
	void save(BS_Out&OUT);
	void set_load(double *F, double *T);
	BS_BallScrew();
	~BS_BallScrew();
};






