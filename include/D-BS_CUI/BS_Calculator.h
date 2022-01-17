#pragma once
#include <stdio.h>
#include <Eigen\Dense>
#include "BS_BallScrew.h"
#include "intel_ode.h"
#include "mkl_rci.h"
#include "Rigid.h"
#include "BS_In.h"
#include "BS_FileIn.h"

class BS_Calculator {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
public:
	//static BS_BallScrew BS;
	static std::shared_ptr<BS_BallScrew> BS;

	// 静解析のために必要なグローバル変数パラメタ群．MKLでは関数の形が決まっているため，このような形で値を受け渡す．
	static struct Stt {
		struct Set {
			int n;				// 入力配列サイズ．
			int m;				// 出力配列サイズ．
			double eps[6];		// ソルバの設定（長さ6の配列）．各要素の意味はマニュアル参照．
			int iter1;			// 最大繰り返し数．
			int iter2;			// 最大試行回数．
			double rs;			// 初期ステップ条件．
			double jac_eps;		// ヤコビアン計算の精度．
		} set[3];
		double v0;		// [m/s]:	シャフト並進速度（//axis）
		double w0;		// [rad/s]:	シャフト回転速度（//axis）
		double wn;		// [rad/s]:	ナット回転速度（//axis）
		int i1;			// ステップ1用．玉番号を格納．
	} stt;

	static struct Dyn {
		struct Set {
			double t_end;		// 計算終了時刻[s]
			double t_step;		// サンプリング時間ステップサイズ（≠最小時間ステップサイズ）
			int nX;				// 玉・内外輪・保持器各要素の状態量（座標・速度・加速度）
			double h;			// 最小ステップサイズ
			double hmin;		// 許容誤差
			double ep;			// しきい値
			double tr;			// 初期ステップサイズ
			double *dpar;
			int ierr;
			int kd[2];
			int ipar[128];
		} set[2];

		map<double, double>wxt;				// 角速度時間変化(first: 時間[s], second:角速度[rad/s])
		map<double, double>vxt;				// 速度時間変化　(first: 時間[s], second:速度[m/s])
		map<double, vector<double>>ft;		// 荷重時間変化　(first: 時間[s], second:荷重[N]/[Nm])

		struct Load {
			//bool is_change;		// 速度条件を切り替えるかどうか（≒変化ステップ2以上）
			struct Param {
				double t;			// [s]:		速度条件が切り替わる時間
				Vector3d F;		// [N]:		tの瞬間の外部荷重
				Vector3d M;		// [Nm]:		tの瞬間の外部荷重
			};
			vector<Param> param;
			int i0;				// 現在の区間
		} lod;
	}dyn;

public:
	static void init_stt(const BS_FileIn::Static & stt, int ballnum);

private:
	static void Stt_Eq0(int * m, int * n, double * x, double * f);
	static void Stt_Eq1(int * m, int * n, double * x, double * f);
	static void Stt_Eq2(int * m, int * n, double * x, double * f);

public:
	static int Stt_solve(int i, double*x);

public:
	static void init_dyn(const BS_FileIn::Dynamic & dyn, int ballnum);
	   
private:
	static void Dyn_Eq0(int * n, double * t, double * y, double * dydt);
	static void Dyn_Eq1(int *n, double *t0, double *y, double *dydt);
	static void Dyn_void(void);
	//static void get_dwdt(double t, double &dwdt);
	static void get_Load(double t, map<double, vector<double>>ft, double *F, double *T);

	//static double get_dvdt(double t, map<double, double> vxt);
	static double get_dvdt(double t, const map<double, double>& vxt);

public:
	static void Dyn_solve(double * x0, double * t, int i);


};
