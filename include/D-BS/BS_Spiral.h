#pragma once
#include <Eigen\Dense>
#include "Spiral.h"
#include "BS_In.h"

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Matrix3d;

class BS_Spiral {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
private:
	Spiral SP;		// 螺旋クラス

public:
	struct Groove {
		Vector2d eta;	// 溝中心の径方向ズレ [m]．軸側+．ナット正，シャフト負．
		double r;		// 溝半径 [m]．
		double r_inv;	// 溝曲率 [1/m]．
		double sigma;	// 表面粗さ [m]
		double R;		// 曲率R [m]．ループ不変量なのでメンバに保存．
	} GV[2];

public:
	// 委譲するだけで速度律速が出るのは悔しいので，なるべく遅くならないようにインライン展開をする．
	inline Vector3d to_eta(const Vector3d&xyz)	{ return this->SP.to_eta2(xyz); };
	inline Vector3d to_xyz(const Vector3d&eta)	{ return this->SP.to_xyz(eta); };
	inline Matrix3d get_xyz2eta(double theta)	{ return this->SP.get_xyz2eta(theta); };
	inline double   get_nd(void)				{ return this->SP.get_nd(); };
	inline double   get_r(void)					{ return this->SP.get_r(); };

	void Nothing(void);
	void init(const BS_In::Cylinder::Spiral & spiral);
	Vector2d get_rho(double cos_alp, int i);
};




