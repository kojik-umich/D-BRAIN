/******************************************************************************* 
!								"B4P_Ring.cpp"	
!													2019/03/22	[Core-T]	楠崎
!	
!	
!*******************************************************************************/
#include "B4P_Ring.h"
#include<iostream>
void B4P_Ring::Nothing(void) {
}

// 接触楕円をトーラスに沿ってスライスするメソッド．(関数内の座標とベクトルはすべて慣性座標系)
void B4P_Ring::calc_slice
	(
	const Vector3d&p,	// in:	[m]:	接触点位置の中心．（慣性座標系）
	double a,			// in:	[m]:	接触楕円長径．
	int msmax,			// in:	[-]:	接触楕円のスライス数．
	int ig,				// in:	[-]:	溝番号．0 or 1;
	double bl_r,		// in:	[m]:	玉半径
	Vector3d*ps		// out:	[m]:	スライスされた各接触点位置．配列サイズ=msmax．（慣性座標系）
	)
{
	// 内外輪中心から接触点へ向かうベクトルを取得
	Vector3d xp = p - this->x;

	// 現在の軸方向の単位ベクトルを算出
	Vector3d ax = this->get_ax();
	
	// 接触点における周方向ベクトル（軸方向とベクトルxpと直交する単位ベクトル）
	Vector3d cs  = (xp.cross(ax)).normalized();
	
	// 接触点における径方向ベクトル（軸方向に垂直）
	Vector3d er = ax.cross(cs);
	
	// 慣性座標系での，点pが乗る円弧の中心O．（=溝R中心と同一と仮定）
	Vector3d O = this->GV[ig].Rx * ax + this->GV[ig].Rr * er + this->x;

	// 楕円をスライスした際の，距離の間隔．
	double da = 2 * a / msmax;

	// 円の計算はコストがかかるため，接触楕円が十分小さいことから2次近似を行う．
	Vector3d op = p - O;
	double op_norm = op.norm();

	if (op_norm < GV[ig].r) {
		double _t = 0;
		std::cout << "err!!!" << op_norm << " " << GV[ig].r << std::endl;
	}
	Vector3d cn = op / op_norm; // 溝断面に垂直なベクトル
	Vector3d ct = cn.cross(cs);	// 接触断面と並行かつcsと垂直な方向ベクトル
	double Rm = 2 * bl_r * this->GV[ig].r / (bl_r + this->GV[ig].r);
	// 各スライスの中心点 pi の座標を算出
	for (int i = 0; i < msmax; i++) {
		double xi_ = i * da - a + 0.5 * da;		// pを基準とした点p_iのctベクトル方向の距離
		double yi_ = - 0.5 * xi_ * xi_ / Rm;	// pを基準とした点p_iのcnベクトル方向の距離
		ps[i] = p + cn * yi_ + ct * xi_;
	}
	return;
}

// パラメータをメンバ変数に代入
void B4P_Ring::init(const B4P_In::Ring&IN, double pcd) {
	this->pcd  = pcd;					// 玉PCD
	for (int i = 0; i < 2; i++) {
		this->GV[i].Rx = IN.Rox[i];			// pcd から見た軸方向変位(-x側)
		this->GV[i].r = IN.R[i];				// 溝径(-x側)
		this->GV[i].r_inv = 1.0 / IN.R[i];	// 溝径の曲率(-x側)
		this->GV[i].Rr = IN.Rod[i] / 2;		// 内外輪中心から溝中心の距離 [m]
		this->GV[i].h = IN.hedge[i];			// 溝肩高さ[m](新規追加)
	}
	this->sigmap   =IN.rms;					// 粗さrms(新規追加)
	this->m			= IN.m;					// 重量
	this->m_inv		= 1.0 / this->m;
	this->E  = IN.E;							// 剛性
	this->nu = IN.por;
	this->I			= Vector3d(IN.Ix, IN.Iyz, IN.Iyz);
	this->I_inv		= Vector3d(1/IN.Ix, 1/IN.Iyz, 1/IN.Iyz);
	return;
}

// 拘束条件用の物体座標系に変換するための回転行列を求める．
// 拘束条件用座標系：x'方向が軸を向いている座標系．ただし，x'方向には回転しない．
Matrix3d B4P_Ring::get_Rth() {
	Vector3d ex_d = this->get_ax();		// 物体座標系x'方向（軸方向）
	Vector3d ex = Vector3d(1, 0, 0);	// 慣性座標系x方向
	// 移動前と移動後の角度を取る
	double angle = acos(ex.dot(ex_d));
	//std::cout << ex_d << std::endl;
	// 回転軸(外積を求めて正規化)
	Vector3d axis = ex.cross(ex_d).normalized();
	// 回転量と回転軸から回転行列を生成
	Matrix3d Rth = Eigen::AngleAxisd(angle, axis).matrix();

	return Rth;
}



// F, T の入力により x, v, q, w 全ての時間微分値を返すメソッド．拘束条件付き．値を全て無次元量にして返す．
// クォータニオンの拘束条件に Baumgarte の方法を適用
void B4P_Ring::get_dydt_(
	const Vector3d&F,	// in:	[N]:	外部荷重．
	const Vector3d&T,	// in:	[Nm]:	外部トルク．
	double*dydt			// out:	[any]:	全微分値．
) {

	// 姿勢拘束（角速度とクォータニオンの両方に拘束をかける）
	Vector4d dqdt = this->get_dqdt_(this->Rx_const.y(), this->Rx_const.z()) * Rigid::t;
	Vector3d dwdt = this->get_dwdt_(T, this->Rx_const.y(), this->Rx_const.z()) * Rigid::t * Rigid::t;

	// 変位拘束
	Vector3d dxdt = this->get_dxdt() / Rigid::l * Rigid::t;
	Vector3d dvdt = this->get_dvdt(F) / Rigid::l * Rigid::t * Rigid::t;
	dvdt = (this->x_const == true).select(Vector3d::Zero(), dvdt);

	dydt[0] = dxdt.x(); dydt[1] = dxdt.y(); dydt[2] = dxdt.z();
	dydt[3] = dvdt.x(); dydt[4] = dvdt.y(); dydt[5] = dvdt.z();
	dydt[6] = dqdt.w(); dydt[7] = dqdt.x(); dydt[8] = dqdt.y(); dydt[9] = dqdt.z();
	dydt[10] = dwdt.x(); dydt[11] = dwdt.y(); dydt[12] = dwdt.z();
	return;
}

// 角度拘束付きクォータニオンの時間微分式．クォータニオンは加算が定義されていないため，ベクトル配列で計算．
Vector4d B4P_Ring::get_dqdt_(bool wy_const, bool wz_const){
	// Quaterniond のコンストラクタは(w, x, y, z)の順番
	Quaterniond w_half(0.0, 0.5*this->w.x(), 0.5*this->w.y(), 0.5*this->w.z());
	Vector4d qw_half = (w_half * q).coeffs();
	Vector4d dqdt;
	double tan_beta = 0;
	double tan_gamma = 0;
	double tau = 0.1;
	

	// 完全拘束
	if (wy_const && wz_const) {
		dqdt = qw_half;
	}
	// wy方向が拘束されている場合，wy方向と軸方向に加速度成分を持たないようにする．
	// (= wy単位ベクトルと軸方向ベクトルの外積の向きのみに加速するようにする)
	else if (wy_const && !wz_const) {
		// Vector4d のコンストラクタは(x, y, z, w)の順番
		Vector4d _Cq(q.z() + tan_beta * q.x(), -q.w() - tan_beta * q.y(),
			q.x() - tan_beta * q.z(), -q.y() + tan_beta * q.w());
		Vector4d Cq = _Cq * 2;
		double C = 2 * (q.x() * q.z() - q.y() * q.w()) +
			(q.x() * q.x() - q.y() * q.y() - q.z() * q.z() + q.w() * q.w()) * tan_beta;
		// 0割りの判定
		if (Cq.norm() < 1e-20) {
			dqdt = qw_half;
		}
		else {
			dqdt = qw_half - Cq * (qw_half.dot(Cq) + 1 / tau * C) / Cq.dot(Cq);
		}
	}
	// wz方向が拘束されている場合，wz方向と軸方向に加速度成分を持たないようにする．
	// (= wy単位ベクトルと軸方向ベクトルの外積の向きのみに加速するようにする)
	else if (!wy_const && wz_const) {
		// Vector4d のコンストラクタは(x, y, z, w)の順番
		Vector4d _Cq(q.y() - tan_gamma * q.x(), q.x() + tan_gamma * q.y(),
			q.w() + tan_gamma * q.z(), q.y() - tan_gamma * q.w());
		Vector4d Cq = _Cq * 2;
		double C = 2 * (q.x() * q.y() + q.z() * q.w()) +
			(q.x() * q.x() - q.y() * q.y() - q.z() * q.z() + q.w() * q.w()) * tan_gamma;
		// 0割りの判定
		if (Cq.norm() < 1e-20) {
			dqdt = qw_half;
		}
		else {
			dqdt = qw_half - Cq * (qw_half.dot(Cq) + 1 / tau * C) / Cq.dot(Cq);
		}
	}
	// 軸方向の回転のみ拘束．
	else {
		dqdt = qw_half;
	}
	return dqdt;
}

// 角速度の時間微分式．軸自転方向と拘束されている成分は変化しないように設定．
Vector3d B4P_Ring::get_dwdt_(const Vector3d&T, bool wy_const, bool wz_const) {
	Vector3d dwdt;

	// 完全拘束
	if (wy_const && wz_const) {
		dwdt = Vector3d::Zero();
	}
	// wy方向が拘束されている場合，wy方向と軸方向に加速度成分を持たないようにする．
	// (= wy単位ベクトルと軸方向ベクトルの外積の向きのみに加速するようにする)
	else if (wy_const && !wz_const) {
		Vector3d ax = this->get_ax();
		Vector3d ey = Vector3d(0, 1, 0);
		Vector3d ew = ax.cross(ey);
		Vector3d _dwdt = this->get_dwdt(T);
		dwdt = _dwdt.dot(ew) * ew;
	}
	// wz方向が拘束されている場合，wz方向と軸方向に加速度成分を持たないようにする．
	// (= wy単位ベクトルと軸方向ベクトルの外積の向きのみに加速するようにする)
	else if (!wy_const && wz_const) {
		Vector3d ax = this->get_ax();
		Vector3d ez = Vector3d(0, 0, 1);
		Vector3d ew = ax.cross(ez);
		Vector3d _dwdt = this->get_dwdt(T);
		dwdt = _dwdt.dot(ew) * ew;
	}
	// 軸方向の回転のみ拘束．
	else {
		Vector3d dwdt_ = this->get_dwdt(T);
		dwdt = this->remove_ax(dwdt_);
	}
	return dwdt;
}

// 静解析　剛性計算：内輪の座標・姿勢を配列 Xstf に格納
void B4P_InnerRing::get_Xstf(double* Xstf) {

	Xstf[0] = this->x.x() / Rigid::l;
	Xstf[1] = this->x.y() / Rigid::l;
	Xstf[2] = this->x.z() / Rigid::l;

	Vector3d ax = this->get_ax();
	Xstf[3] = ax.y();
	Xstf[4] = ax.z();

	return;
}

// 静解析　剛性計算：内輪の座標，姿勢を入力
void B4P_InnerRing::set_Xstf(double* Xstf) {

	double x0 = sqrt(1.0 - Numeric::Square(Xstf[3]) - Numeric::Square(Xstf[4]));
	Vector3d ax(x0, Xstf[3], Xstf[4]);
	this->set_ax(ax);
	Vector3d ax_ = this->get_ax();

	this->x = Vector3d(Xstf[0], Xstf[1], Xstf[2]) * Rigid::l;

	return;
}

// 位相角から慣性座標系のベクトルを溝直行断面座標系へ変換するメソッド．
Vector3d B4P_Ring::ine_to_XZvector
(							// out:	[any]:	0成分：リング位相角．1成分：X座標．2成分：Z座標．
	const Vector3d&a,		// in : [any]:	変換したいベクトル．（慣性座標系）
	const Matrix3d&xyz2XYZ	// in : [rad]:	位相角分の回転行列（リング座標系）
) {
	Vector3d xyz = this->to_myvector(a);
	Vector3d XYZ = xyz2XYZ * xyz;
	return XYZ;
}

// リング座標系から溝直行座標系へ変換する行列を取得するメソッド．
Matrix3d B4P_Ring::get_xyz2XYZ
(						// out:	[any]:	0成分：リング位相角．1成分：X座標．2成分：Z座標．
	double th			// in : [rad]:	位相角（リング座標系）
) {
	double cos_th = cos(th);
	double sin_th = sin(th);
	Matrix3d xyz2XYZ;
	xyz2XYZ <<
		1.0, 0.0, 0.0,
		0.0, cos_th, sin_th,
		0.0, -sin_th, cos_th;
	return xyz2XYZ;
}

// 慣性座標系からリングの溝直交座標系に変換するメソッド．
Vector3d B4P_Ring::ine_to_XZcoord
(						// out:	[rad],[m]:	0成分：リング位相角．1成分：X座標．2成分：Z座標．
	const Vector3d& x	// in : [m]:		慣性座標系の変換したい座標．
)
{
	Vector3d xyz = this->to_mycoord(x);
	double th = atan2(-xyz[1], xyz[2]);
	double X = xyz[0];
	double Z = sqrt(xyz[1] * xyz[1] + xyz[2] * xyz[2]) - 0.5 * this->pcd;
	Vector3d thXZ = Vector3d(th, X, Z);
	return thXZ;
}

// リングの溝直交座標系から慣性座標系に変換するメソッド．
Vector3d B4P_Ring::XZ_to_inecoord
(						// out:	[m]:		慣性座標系の座標．
	const Vector3d& thXZ	// in :	[rad],[m]:	変換したいリングの溝直交座標系．
)
{
	double th = thXZ[0];
	double x = thXZ[1];
	double R = 0.5 * this->pcd + thXZ[2];
	double y = -R * sin(th);
	double z = R * cos(th);
	Vector3d xyz = Vector3d(x, y, z);
	Vector3d xyz_ = this->to_inecoord(xyz);
	return xyz_;
}

// 与えられた任意のベクトルを，ある位相角での溝直角断面上の2次元ベクトルに射影する．
Vector2d B4P_Ring::to_cross_sction
(					// out:	[Any]:	入力と同じ次元．第0成分：軸方向，第1成分：径方向．
	const Vector3d&v,	// in:	[Any]:	任意のベクトル（慣性座標系）．
	double th			// in:	[rad]:	中立位置で見たときのリング位相角．
)
{
	Vector3d v_ = this->to_myvector(v);
	return Vector2d(v_[0], -sin(th)*v_[1] + cos(th)*v_[2]);
}


// 常微分方程式の解から物体の各パラメータを代入．
// クォータニオンは数値誤差があるため，数値代入時にも拘束条件を適用．
void B4P_Ring::set_y_(const Vector3d&x, const Vector3d&v, const Quaterniond&q, const Vector3d&w) {
	
	this->x = x * Rigid::l;
	this->v = v * Rigid::l / Rigid::t;
	this->w = w / Rigid::t;
	Quaterniond _q = q;
	// クォータニオンは数値の計算誤差があるため，軸方向ベクトルが拘束されている向きに動かないようにする

	if (this->Rx_const[1] == true && _q.w() != 0) {
		_q.y() = _q.x() * _q.z() / _q.w();
	}
	if (this->Rx_const[2] == true && _q.w() != 0) {
		_q.z() = -_q.x() * _q.y() / _q.w();
	}
	this->q = _q;

	this->q.normalize();

}









//
//Vector3d B4P_Ring::to_myXYZ
//	(					// out:	[Any]:	入力と同じ次元．第0成分：軸方向，第1成分：径方向．
//	const Vector3d&a	// in:	[Any]:	任意のベクトル（慣性座標系）．
//	)
//{
//	Vector3d u_  = Vector3d(0.0, -sin(th), cos(th));	// 外輪座標系での進行方向ベクトル．
//	Vector3d u__ = this->to_inevelocity(u_);			// 慣性座標系への変換．
//	return u__;
//}
//
//void B4P_Ring::init(const B4P_In&FI, const InnerOuter&innerouter){
//	switch (innerouter)
//	{
//	case B4P_Ring::Inner:
//		this->pcd  = (FI.igrv_r_trs[0] + FI.igrv_r_trs[1]) / 2 * 2;
//		this->GV[0].x = -FI.X_disti/2;
//		this->GV[1].x = FI.X_disti/2;
//		this->GV[0].r = FI.rgi[0];
//		this->GV[1].r = FI.rgi[1];
//		this->GV[0].R = FI.igrv_r_trs[0];
//		this->GV[1].R = FI.igrv_r_trs[1];
//		this->m  = Numeric::pi/6*pcd*pcd*pcd*FI.rhop[1];
//		this->E  = FI.Ep[1];
//		this->nu = FI.etap[1];
//		this->m  = Numeric::pi/6*pcd*pcd*pcd*FI.rhop[1];	// リングだけど球体で仮定．
//
//	case B4P_Ring::Outer:
//		this->pcd  = (FI.ogrv_r_trs[0] + FI.ogrv_r_trs[1]) / 2 * 2;
//		this->GV[0].x = -FI.X_disto/2;
//		this->GV[1].x = FI.X_disto/2;
//		this->GV[0].r = FI.rgo[0];
//		this->GV[1].r = FI.rgo[1];
//		this->GV[0].R = FI.ogrv_r_trs[0];
//		this->GV[1].R = FI.ogrv_r_trs[1];
//		this->m  = Numeric::pi/6*pcd*pcd*pcd*FI.rhop[2];
//		this->E  = FI.Ep[2];
//		this->nu = FI.etap[2];
//		this->m  = Numeric::pi/6*pcd*pcd*pcd*FI.rhop[2];	// リングだけど球体で仮定．
//
//	default:
//		break;
//	}
//
//	double I_= this->m*pcd*pcd/10;
//	this->I << I_, I_, I_;
//	this->m_inv = 1.0 / this->m;
//	this->I_inv << 1 / I_, 1 / I_, 1 / I_;
//}
//$$Outer $$Inner[内外輪],0密度,1ヤング率,2por,3粗さrms,4溝R(-x),5溝R(+x),6中心O(-x),7中心O(+x),8溝R中心直径(-x),9溝R中心直径(+x),10肩高さ(-x), 11肩高さ(+x), 12玉PCD
//void B4P_Ring::init(const VectorXd &IOringparam, double pcd){
//	//switch (innerouter)
//	//{
//	//case B4P_Ring::Inner:
//		//this->pcd  = IOringparam[12]; // 玉PCD
//		this->pcd  = pcd; // 玉PCD
//		this->GV[0].x = IOringparam[6];	// pcd から見た軸方向変位(-x側)
//		this->GV[1].x = IOringparam[7];	// pcd から見た軸方向変位(+x側)
//		this->GV[0].r = IOringparam[4];	// 溝径(-x側)
//		this->GV[1].r = IOringparam[5];	// 溝径(+x側)
//		this->GV[0].R = IOringparam[8];	// 溝中心pcd [m]
//		this->GV[1].R = IOringparam[9];	// 溝中心pcd [m]
//		this->GV[0].h = IOringparam[10];	// 溝肩高さ[m](新規追加)
//		this->GV[1].h = IOringparam[11];	// 溝肩高さ[m](新規追加)
//		this->sigmap   = IOringparam[3];	// 粗さrms(新規追加)
//		//this->m  = Numeric::pi/6*pcd*pcd*pcd*FI.rhop[1];//重量自動計算
//		this->E  = IOringparam[1];
//		this->nu = IOringparam[2];
//		this->m  = Numeric::pi/6*pcd*pcd*pcd*IOringparam[0];	// リングだけど球体で仮定????
//
//	//case B4P_Ring::Outer:
//		//this->pcd  = (IOringparam[8] + IOringparam[9]) / 2 * 2;
//
//		//this->GV[0].x = -(IOringparam[6]+IOringparam[7])/2;
//		//this->GV[1].x = (IOringparam[6]+IOringparam[7])/2;
//		//this->GV[0].r = IOringparam[4];
//		//this->GV[1].r = IOringparam[5];
//		//this->GV[0].R = IOringparam[8];
//		//this->GV[1].R = IOringparam[9];
//		//this->GV[0].h = IOringparam[10];
//		//this->GV[1].h = IOringparam[11];
//		////??  = IOringparam[3];//粗さrms
//		////this->m  = Numeric::pi/6*pcd*pcd*pcd*FI.rhop[1];//重量自動計算
//		//this->E  = IOringparam[1];
//		//this->nu = IOringparam[2];
//
//		//this->m  = Numeric::pi/6*pcd*pcd*pcd*IOringparam[0];	// リングだけど球体で仮定．
//
//	//default:
//	//	break;
//	//}
//
//	double I_= this->m*pcd*pcd/10;
//	this->I << I_, I_, I_;
//	this->m_inv = 1.0 / this->m;
//	this->I_inv << 1 / I_, 1 / I_, 1 / I_;
//}


//// 慣性座標系速度からリングの溝直交座標系速度へ変換するメソッド．
//double B4P_Ring::ine_to_Yvelocity
//	(					// out:	[m/s]:	リングの溝直交座標系における速度1成分．
//	const Vector3d& v,	// in : [m/s]:	慣性座標系の変換したい速度．
//	double th			// in:	[rad]:	リング位相角．
//	)
//{
//	Vector3d u = this->get_cd(th);	// 位相角での進行方向ベクトル．（慣性座標系）
//	double   v_= u.dot(v);					// 速度の進行方向成分．
//	return   v_;
//}
//
//// 慣性座標系速度からリングの溝直交座標系速度へ変換するメソッド．
//Vector3d B4P_Ring::Y_to_inevelocity
//	(				// out:	[m/s]:	慣性座標系における速度3成分．
//	double v,		// in : [m/s]:	リングの溝直交座標系の変換したい速度．
//	double th		// in:	[rad]:	リング位相角．
//	)
//{
//	Vector3d u  = this->get_cd(th);
//	Vector3d v_ = u * v;
//	return v_;
//}
//

//
//
//// ある位相角での周方向ベクトル(CircumferentialDirection)を取得するメソッド．（q=[1,0,0,0]， θ = 45° で [0,-1/√2,-1/√2] となる）
//Vector3d B4P_Ring::get_cd
//	(				// out:	[m/s]:	慣性座標系における速度3成分．
//	double th		// in:	[rad]:	中立位置で見たときのリング位相角．
//	)
//{
//	Vector3d u_  = Vector3d(0.0, -cos(th), -sin(th));	// リング座標系での進行方向ベクトル．
//	Vector3d u__ = this->to_inevelocity(u_);			// 慣性座標系への変換．
//	return u__;
//}
//
//// ある位相角での径方向ベクトル(RadialDirection)を取得するメソッド．（q=[1,0,0,0]， θ = 45° で [0,-1/√2,1/√2] となる）
//Vector3d B4P_Ring::get_rd
//	(				// out:	[m/s]:	慣性座標系における速度3成分．
//	double th		// in:	[rad]:	中立位置で見たときのリング位相角．
//	)
//{
//	Vector3d u_  = Vector3d(0.0, -sin(th), cos(th));	// 外輪座標系での進行方向ベクトル．
//	Vector3d u__ = this->to_inevelocity(u_);			// 慣性座標系への変換．
//	return u__;
//}


