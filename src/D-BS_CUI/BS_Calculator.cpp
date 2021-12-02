/*******************************************************************************
!								"BS_Calculator.cpp"
!													2020/04/10	[Core-T]	楠崎
!
! このクラスを使用するにはMKL NonLinear Solverのdtrnlsp_solveと
! MKL ODE Solverのodeintにパスを通す必要があります．
! マニュアルの「MKLのインストールの方法」を終えたのち，以下の設定を行って下さい．
!
! 静解析編
! 1. インクルードディレクトリに"$(ICPP_COMPILER18)mkl\include"を追加．
! 2. ライブラリディレクトリに"$(ICLibDir)"と"$(ICInstallDir)mkl\lib\ia32_win"を追加．
! 3. リンカーの入力に"mkl_sequential.lib"と"mkl_intel_c.lib"と"mkl_rt.lib"と"mkl_core.lib"を追加．
!
! 動解析編
! 1. ライブラリディレクトリに"$(ProjectDir)intel_ode\lib\ia32"を追加．
! 2. リンカーの入力に"libiode_ia32.lib"を追加．
!
!*******************************************************************************/

#include "BS_Calculator.h"

// このオブジェクト内でのグローバル変数的役割．
BS_BallScrew BS_Calculator::BS = BS_BallScrew();
BS_Calculator::Stt BS_Calculator::stt;
BS_Calculator::Dyn BS_Calculator::dyn;

const int shaft_DOF = 5;

// 静解析の設定初期化．
void BS_Calculator::init_stt(const BS_FileIn::Static&stt, int ballnum) {

	// 静解析設定
	for (int i = 0; i < 3; i++) {
		// ソルバの設定（長さ6の配列）．各要素の意味はマニュアル参照．
		BS_Calculator::stt.set[i].eps[0] = stt.set[i].eps[0];
		BS_Calculator::stt.set[i].eps[1] = Numeric::Square(stt.set[i].eps[1]);
		BS_Calculator::stt.set[i].eps[2] = Numeric::Square(stt.set[i].eps[2]);
		BS_Calculator::stt.set[i].eps[3] = Numeric::Square(stt.set[i].eps[3]);
		BS_Calculator::stt.set[i].eps[4] = Numeric::Square(stt.set[i].eps[4]);
		BS_Calculator::stt.set[i].eps[5] = stt.set[i].eps[5];

		// 最大繰り返し数．
		BS_Calculator::stt.set[i].iter1 = stt.set[i].iter1;

		// 最大試行回数．
		BS_Calculator::stt.set[i].iter2 = stt.set[i].iter2;

		// 初期ステップ条件．
		BS_Calculator::stt.set[i].rs = stt.set[i].rs;

		// ヤコビアン計算の精度．
		BS_Calculator::stt.set[i].jac_eps = stt.set[i].jac_eps;
	}
	// 入力配列サイズ．
	BS_Calculator::stt.set[0].n = ballnum * 2 + shaft_DOF;
	BS_Calculator::stt.set[1].n = ballnum * 2 + shaft_DOF;
	BS_Calculator::stt.set[2].n = ballnum * 6 + shaft_DOF;

	// 出力配列サイズ．= 入力サイズ
	for (int i = 0; i < 3; i++)
		BS_Calculator::stt.set[i].m = BS_Calculator::stt.set[i].n;

	BS_Calculator::stt.v0 = stt.v0;
	BS_Calculator::stt.w0 = stt.w0;
	BS_Calculator::stt.wn = stt.wn;

	return;
}

// MKL "dtrnlsp_solve" を使った静解析ソルバ．使い方はマニュアル参照．（変数のアサインをマニュアルと統一させています．）
int BS_Calculator::Stt_solve
(
	int i,			// in:		[-]:	静解析のstep番号．Stt_Eq*が各iに対応．
	double*x		// inout:	[-]:	解きたい関数の引数．
) {
	// ソルバの初期化．
	_TRNSP_HANDLE_t handle;

	// ここは後で設定に反映．
	double*x_lw = new double[BS_Calculator::stt.set[i].n];
	double*x_up = new double[BS_Calculator::stt.set[i].n];

	for (int j = 0; j < BS_Calculator::stt.set[i].n; j++) {
		x_lw[j] = -1e-4;
		x_up[j] = 1e-4;
	}

	//// ちょっと試す
	////double*f0 = new double[BS_Calculator::stt.set[i].n];
	//double*x0 = new double[BS_Calculator::stt.set[i].n];
	//double*dx0 = new double[BS_Calculator::stt.set[i].n];

	//for (int j = 0; j < BS_Calculator::stt.set[i].n; j++) {
	//	//f0[j] = 0.0;
	//	x0[j] = 0.0;
	//	dx0[j] = 0.0;
	//}

	dtrnlspbc_init(&handle, &stt.set[i].n, &stt.set[i].m, x, x_lw, x_up, stt.set[i].eps, &stt.set[i].iter1, &stt.set[i].iter2, &stt.set[i].rs);

	// 収束計算のループ．関数値およびヤコビアン配列の動的確保．
	double *fvec = new double[stt.set[i].m];
	double *fjac = new double[stt.set[i].m * stt.set[i].n];
	int RCI_Request = 0;

	while (true) {

		dtrnlspbc_solve(&handle, fvec, fjac, &RCI_Request);

		//for (int j = 0; j < BS_Calculator::stt.set[i].n; j++) {
		//	double G = x[j] - x0[j];
		//	dx0[j] = 0.5 * dx0[j] + 0.3 * G;
		//	x[j] = x0[j] + dx0[j];
		//	x0[j] = x[j];
		//}

		if (RCI_Request == -1) {
			printf("最大繰り返し数 %d を超過しました．収束に失敗したためプログラムを終了します．\n", stt.set[0].iter1);
			break;
		}
		else if (RCI_Request == -2) {
			printf("修正量が %e [μm] を下回りました．収束に失敗したためプログラムを終了します．\n", stt.set[i].eps[0]);
			break;
		}
		else if (RCI_Request == -3) {
			printf("部材の力の釣り合いが基準値 %e [N] を下回りました．計算は収束しました．\n", sqrt(stt.set[i].eps[1]));
			break;
		}
		else if (RCI_Request == -4) {
			printf("ヤコビアンの行列式が %e [N] を下回りました．収束に失敗したためプログラムを終了します．\n", sqrt(stt.set[i].eps[2]));
			break;
		}
		else if (RCI_Request == -5) {
			printf("探索域が %e [μm] を下回りました．収束に失敗したためプログラムを終了します．\n", sqrt(stt.set[i].eps[3]));
			break;
		}
		else if (RCI_Request == -6) {
			printf("1ステップの変化量が %e [N] を下回りました．収束に失敗したためプログラムを終了します．\n", sqrt(stt.set[i].eps[4]));
			break;
		}
		else if (RCI_Request == 1) {
			if (i == 0) {
				Stt_Eq0(&stt.set[i].m, &stt.set[i].n, x, fvec);
			}
			else if (i == 1) {
				Stt_Eq1(&stt.set[i].m, &stt.set[i].n, x, fvec);
			}
			else if (i == 2) {
				Stt_Eq2(&stt.set[i].m, &stt.set[i].n, x, fvec);
			}
		}
		else if (RCI_Request == 2) {
			if (i == 0) {
				djacobi(Stt_Eq0, &stt.set[i].n, &stt.set[i].m, fjac, x, &stt.set[i].jac_eps);
			}
			else if (i == 1) {
				djacobi(Stt_Eq1, &stt.set[i].n, &stt.set[i].m, fjac, x, &stt.set[i].jac_eps);
			}
			else if (i == 2) {
				djacobi(Stt_Eq2, &stt.set[i].n, &stt.set[i].m, fjac, x, &stt.set[i].jac_eps);
			}
		}
	}
	// ソルバのシャットダウン．
	dtrnlspbc_delete(&handle);

	// 配列の開放．
	delete[] fvec, fjac;
	MKL_Free_Buffers();

	return RCI_Request;
}

// 静解析のSTEP0，剛性計算．変位を入力，荷重を出力に見立ててボールねじをブラックボックスとして扱う．
void BS_Calculator::Stt_Eq0(
	int    *m, 		// in:	[-]:	出力配列サイズ．
	int    *n, 		// in:	[-]:	入力配列サイズ．
	double *x,		// in:	[-]:	関数の引数．
	double *f		// out:	[-]:	関数の戻り値．
) {
	BS_Calculator::BS.set_y0(x, stt.v0, stt.w0);
	BS_Calculator::BS.get_F0(f);
	return;
}

// 静解析のSTEP1，滑り計算．各玉ごとに，純転がりする値を求める．
void BS_Calculator::Stt_Eq1(
	int    *m, 		// in:	[-]:	出力配列サイズ．
	int    *n, 		// in:	[-]:	入力配列サイズ．
	double *x,		// in:	[-]:	関数の引数．
	double *f		// out:	[-]:	関数の戻り値．
) {
	BS_Calculator::BS.set_y0(x, stt.v0, stt.w0);
	BS_Calculator::BS.get_F1(f);
	return;
}

//
// 静解析のSTEP2，自公転釣り合い計算．
void BS_Calculator::Stt_Eq2(
	int    *m, 		// in:	[-]:	出力配列サイズ．
	int    *n, 		// in:	[-]:	入力配列サイズ．
	double *x,		// in:	[-]:	関数の引数．
	double *f		// out:	[-]:	関数の戻り値．
) {
	BS_Calculator::BS.set_y2(x);
	BS_Calculator::BS.get_F2(f);
	return;
}

// 動解析の設定初期化．
void BS_Calculator::init_dyn(const BS_FileIn::Dynamic&dyn, int ballnum) {

	// 動解析設定
	BS_Calculator::dyn.set[0].nX = 4 * ballnum + 10;
	BS_Calculator::dyn.set[1].nX = 13 * (ballnum + 2);

	for (int i = 0; i < 2; i++) {

		int num_loop = dyn.set[i].sampling;

		BS_Calculator::dyn.set[i].t_end = dyn.set[i].calctime;

		BS_Calculator::dyn.set[i].t_step = dyn.set[i].calctime / num_loop;
		BS_Calculator::dyn.set[i].h = dyn.set[i].h;
		BS_Calculator::dyn.set[i].hmin = dyn.set[i].hmin;
		BS_Calculator::dyn.set[i].ep = dyn.set[i].ep;
		BS_Calculator::dyn.set[i].tr = dyn.set[i].tr;
		BS_Calculator::dyn.set[i].dpar = new double[13 * BS_Calculator::dyn.set[i].nX];

		for (int j = 0; j < 128; j++)
			BS_Calculator::dyn.set[i].ipar[j] = 0;

		BS_Calculator::dyn.set[i].ipar[1] = 1;
	}

	BS_Calculator::dyn.wxt = dyn.wxt;
	BS_Calculator::dyn.vxt = dyn.vxt;
	BS_Calculator::dyn.ft = dyn.ft;

	return;
}

// ベクトル表現による運動方程式の定義．12*2+9*n次元で構成されている．
// intel_odeライブラリの書式に準じているため，詳細はそちらのマニュアル参照．
void BS_Calculator::Dyn_Eq0(int*n, double*t, double*y, double*dydt) {

	BS.set_dyn_y0(y);
	BS.get_dyn_dydt0(dydt);

	//for (size_t i = 0; i < 5; i++) 
	//	std::cout << "\t" << y[i];
	//std::cout << std::endl;

	return;
}

// ベクトル表現による運動方程式の定義．12*2+9*n次元で構成されている．
// intel_odeライブラリの書式に準じているため，詳細はそちらのマニュアル参照．
void BS_Calculator::Dyn_Eq1(int*n, double*t, double*y, double*dydt) {
	// 時間から入力条件を取得
	double t0 = *t * Rigid::t;
	double dvdt = 0.0;
	double dwdt = 0.0;
	double F[3], T[3];
	BS_Calculator::get_Load(t0, dyn.ft, F, T);
	//BS_Calculator::get_dwdt(t0, dwdt);
	dwdt = BS_Calculator::get_dvdt(t0, dyn.wxt);
	dvdt = BS_Calculator::get_dvdt(t0, dyn.vxt);

	// まず入力値から各部材の変位を更新．
	BS.set_dyn_y1(y);
	BS.set_load(F, T);
	// 更新した変位から荷重を計算．
	BS.get_dyn_dydt1(dydt, dvdt, dwdt);

	return;
}

// ダミー関数．引数に必要だが使わないため適当に定義した．
void BS_Calculator::Dyn_void(void) {
	return;
}

// Intel ODE Solver を用いて動解析を行う．（内部の記載はほとんどコピペ．dodesolは取得後変更無し．）
void BS_Calculator::Dyn_solve(double*x0, double*t, int i) {

	double t0 = t[0] / Rigid::t;
	double t1 = t[1] / Rigid::t;

	if (i == 0)
		dodesol(dyn.set[i].ipar, &dyn.set[i].nX, &t0, &t1, x0, &Dyn_Eq0, &Dyn_void, &dyn.set[i].h, &dyn.set[i].hmin, &dyn.set[i].ep, &dyn.set[i].tr, dyn.set[i].dpar, dyn.set[i].kd, &dyn.set[i].ierr);

	else if (i == 1)
		dodesol(dyn.set[i].ipar, &dyn.set[i].nX, &t0, &t1, x0, &Dyn_Eq1, &Dyn_void, &dyn.set[i].h, &dyn.set[i].hmin, &dyn.set[i].ep, &dyn.set[i].tr, dyn.set[i].dpar, dyn.set[i].kd, &dyn.set[i].ierr);


	t[0] = t0 * Rigid::t;
	t[1] = t1 * Rigid::t;
	
	return;
}

// w・vの時間変化から，時刻tにおける加速度を計算
double BS_Calculator::get_dvdt(double t, const map<double, double> &vxt) {

	// tが含まれる区間の上限を導出
	auto upper = vxt.upper_bound(t);

	// 区間の上限をはみ出した場合と，時間変化ステップ数が1以下の時は時間変化なし
	if (upper == vxt.end() || vxt.size() <= 1) {
		return 0;
	}
	double t1 = upper->first;
	double vd1 = upper->second;

	--upper;
	double t0 = upper->first;
	double vd0 = upper->second;

	double axd = (vd1 - vd0) / (t1 - t0);
	return axd;

}


// 時刻tにおける外部荷重を線形補間で計算
void BS_Calculator::get_Load(double t, map<double, vector<double>>ft, double *F, double *T) {

	// 荷重条件が入力されていない場合は荷重入力しない
	if (ft.size() == 0) {
		auto min = ft.begin();
		for (int i = 0; i < 3; i++) {
			F[i] = 0;
			T[i] = 0;
		}
		return;
	}

	// tが含まれる区間の上限を導出
	auto upper = ft.upper_bound(t);

	// 区間の上限をはみ出した場合はtが最大の数値を取得
	if (upper == ft.end()) {
		auto max = ft.rbegin();
		for (int i = 0; i < 3; i++) {
			F[i] = max->second[i];
			T[i] = max->second[i + 3];
		}
		return;
	}

	// 線形補間
	double t1 = upper->first;
	vector<double> FT1 = upper->second;
	vector<double> F1{ FT1.begin(), FT1.begin() + 3 };
	vector<double> T1{ FT1.begin() + 3, FT1.begin() + 6 };
	--upper;

	double t0 = upper->first;
	vector<double> FT0 = upper->second;
	vector<double> F0{ FT0.begin(), FT0.begin() + 3 };
	vector<double> T0{ FT0.begin() + 3, FT0.begin() + 6 };
	for (int i = 0; i < 3; i++) {
		F[i] = F0[i] + (F1[i] - F0[i]) * (t - t0) / (t1 - t0);
		T[i] = T0[i] + (T1[i] - T0[i]) * (t - t0) / (t1 - t0);
	}
	return;
}


/*
int main(void) {
	// 入力配列サイズ．
	int n = 4;

	// 出力配列サイズ．
	int m = 4;

	// ソルバの設定（長さ6の配列）．各要素の意味はマニュアル参照．
	double eps[6] ={ 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5 };

	// 最大繰り返し数．
	int iter1 = 1e3;

	// 最大試行回数．
	int iter2 = 1e2;

	// 初期ステップ条件．
	double rs = 0.0;

	// ヤコビアン計算の精度．
	double jac_eps = 1.0e-8;

	// 関数の引数（の初期値設定も含む．）
	double *x = new double[n];
	init_x(x, n);

	// MKL "dtrnlsp_solve" を使って解く．
	Stt_solve(x, n, m, eps, iter1, iter2, rs, jac_eps);

	delete[] x;

	return 0;
}
*/
//void BS_Calculator::initp(const BS_In&IN) {
//
//	// ソルバの設定（長さ6の配列）．各要素の意味はマニュアル参照．
//	for (int i = 0; i < 6; i++)
//		sttsetp.eps[i] = 1e-30;
//	sttsetp.eps[1] = Numeric::Square(IN.preload.F_lim);
//
//	// 最大繰り返し数．
//	sttsetp.iter1 = 10000;
//
//	// 最大試行回数．
//	sttsetp.iter2 = 100;
//
//	// 初期ステップ条件．
//	sttsetp.rs = 0.0;
//
//	// ヤコビアン計算の精度．
//	sttsetp.jac_eps = 1e-8;
//
//	// 入力配列サイズ．
//	sttsetp.n = IN.ballnum * 2 + 2;
//
//	// 出力配列サイズ．= 入力サイズ
//	sttsetp.m = sttsetp.n;
//
//	return;
//}
//
//// 静解析のSTEP0，剛性計算．変位を入力，荷重を出力に見立ててボールねじをブラックボックスとして扱う．
//void BS_Calculator::Stt_Eqp(
//	int    *m, 		// in:	[-]:	出力配列サイズ．
//	int    *n, 		// in:	[-]:	入力配列サイズ．
//	double *x,		// in:	[-]:	関数の引数．
//	double *f		// out:	[-]:	関数の戻り値．
//) {
//	BS_Calculator::BS.set_yp(x);
//	BS_Calculator::BS.get_Fp(f);
//	return;
//}
//
//// MKL "dtrnlsp_solve" を使った静解析ソルバ．使い方はマニュアル参照．（変数のアサインをマニュアルと統一させています．）
//void BS_Calculator::Stt_solvep
//(
//	double*x			// inout:	[-]:	解きたい関数の引数．
//) {
//	int n          = sttsetp.n;			// [-]:	入力配列サイズ．
//	int m          = sttsetp.m;			// [-]:	出力配列サイズ．
//	double eps[6];
//	for (int i = 0; i < 6; i++)
//		eps[i]     = sttsetp.eps[i];	// [-]:	ソルバの設定（長さ6の配列）．各要素の意味はマニュアル参照
//	int iter1      = sttsetp.iter1;		// [-]:	最大繰り返し数．
//	int iter2      = sttsetp.iter2;		// [-]:	最大試行回数．
//	double rs      = sttsetp.rs;		// [-]:	初期ステップ条件．
//	double jac_eps = sttsetp.jac_eps;	// [-]:	ヤコビアン計算の精度．
//
//	// ソルバの初期化．
//	_TRNSP_HANDLE_t handle;
//	dtrnlsp_init(&handle, &n, &m, x, eps, &iter1, &iter2, &rs);
//
//	// 収束計算のループ．関数値およびヤコビアン配列の動的確保．
//	double *fvec = new double[m];
//	double *fjac = new double[m * n];
//	int RCI_Request = 0;
//
//	while (true) {
//		dtrnlsp_solve(&handle, fvec, fjac, &RCI_Request);
//
//		if (RCI_Request == -1) {
//			printf("最大繰り返し数 %d を超過しました．収束に失敗しました．\n", iter1);
//			break;
//		}
//		if (RCI_Request == -2) {
//			printf("修正量が %e [μm] を下回りました．極小値への収束です．\n", eps[0]);
//			break;
//		}
//		if (RCI_Request == -3) {
//			printf("部材の力の釣り合いが基準値 %e [N] を下回りました．計算は収束しました．\n", sqrt(eps[1]));
//			break;
//		}
//		if (RCI_Request == -4) {
//			printf("ヤコビアンの行列式が %e [N] を下回りました．収束に失敗しました．\n", sqrt(eps[2]));
//			break;
//		}
//		if (RCI_Request == -5) {
//			printf("探索域が %e [μm] を下回りました．極小値への収束です．\n", sqrt(eps[3]));
//			break;
//		}
//		if (RCI_Request == -6) {
//			printf("1ステップの変化量が %e [N] を下回りました．極小値への収束です．\n", sqrt(eps[4]));
//			break;
//		}
//		if (RCI_Request == 1)
//			BS_Calculator::Stt_Eqp(&m, &n, x, fvec);
//
//		if (RCI_Request == 2)
//			djacobi(BS_Calculator::Stt_Eqp, &n, &m, fjac, x, &jac_eps);
//	}
//
//	// ソルバのシャットダウン．
//	dtrnlsp_delete(&handle);
//
//	// 配列の開放．
//	delete[] fvec, fjac;
//	MKL_Free_Buffers();
//
//	return;
//}
	//// 与圧の計算設定を適当に設定．
	//BS_Calculator::initp(IN);


//// 静解析のSTEP1，滑り計算．各玉ごとに，純転がりする値を求める．
//void BS_Calculator::Stt_Eq2(
//	int    *m, 		// in:	[-]:	出力配列サイズ．
//	int    *n, 		// in:	[-]:	入力配列サイズ．
//	double *x,		// in:	[-]:	関数の引数．
//	double *f		// out:	[-]:	関数の戻り値．
//) {
//	BS_Calculator::BS.set_y2(x);
//	BS_Calculator::BS.get_F2(f);
//	return;
//}
//
//// MKL "dtrnlsp_solve" を使った静解析ソルバ．使い方はマニュアル参照．（変数のアサインをマニュアルと統一させています．）
//void BS_Calculator::Stt_solve2
//(
//	double*x			// inout:	[-]:	解きたい関数の引数．
//) {
//	int n          = stt.set[2].n;			// [-]:	入力配列サイズ．
//	int m          = stt.set[2].m;			// [-]:	出力配列サイズ．
//	double eps[6];
//	for (int i = 0; i < 6; i++)
//		eps[i]     = stt.set[2].eps[i];		// [-]:	ソルバの設定（長さ6の配列）．各要素の意味はマニュアル参照
//	int iter1      = stt.set[2].iter1;		// [-]:	最大繰り返し数．
//	int iter2      = stt.set[2].iter2;		// [-]:	最大試行回数．
//	double rs      = stt.set[2].rs;			// [-]:	初期ステップ条件．
//	double jac_eps = stt.set[2].jac_eps;		// [-]:	ヤコビアン計算の精度．
//
//	// ソルバの初期化．
//	_TRNSP_HANDLE_t handle;
//	dtrnlsp_init(&handle, &n, &m, x, eps, &iter1, &iter2, &rs);
//
//	// 収束計算のループ．関数値およびヤコビアン配列の動的確保．
//	double *fvec = new double[m];
//	double *fjac = new double[m * n];
//	int RCI_Request = 0;
//
//	while (true) {
//		dtrnlsp_solve(&handle, fvec, fjac, &RCI_Request);
//
//		if (RCI_Request == -1) {
//			printf("最大繰り返し数 %d を超過しました．収束に失敗しました．\n", iter1);
//			break;
//		}
//		if (RCI_Request == -2) {
//			printf("修正量が %e [m/s] を下回りました．玉に滑りが発生しています．\n", eps[0] * Rigid::l / Rigid::t);
//			break;
//		}
//		if (RCI_Request == -3) {
//			printf("部材の力の釣り合いが基準値 %e [N] を下回りました．計算は収束しました．\n", sqrt(eps[1]));
//			break;
//		}
//		if (RCI_Request == -4) {
//			printf("ヤコビアンの行列式が %e [N] を下回りました．収束に失敗しました．\n", sqrt(eps[2]));
//			break;
//		}
//		if (RCI_Request == -5) {
//			printf("探索域が %e [m/s] を下回りました．玉に滑りが発生しています．\n", sqrt(eps[3])* Rigid::l / Rigid::t);
//			break;
//		}
//		if (RCI_Request == -6) {
//			printf("1ステップの変化量が %e [N] を下回りました．玉に滑りが発生しています．\n", sqrt(eps[4]));
//			break;
//		}
//		if (RCI_Request == 1)
//			BS_Calculator::Stt_Eq2(&m, &n, x, fvec);
//
//		if (RCI_Request == 2)
//			djacobi(BS_Calculator::Stt_Eq2, &n, &m, fjac, x, &jac_eps);
//	}
//
//	// ソルバのシャットダウン．
//	dtrnlsp_delete(&handle);
//
//	// 配列の開放．
//	delete[] fvec, fjac;
//	MKL_Free_Buffers();
//
//	return;
//}



// MKL "dtrnlsp_solve" を使った静解析ソルバ．使い方はマニュアル参照．（変数のアサインをマニュアルと統一させています．）
//int BS_Calculator::Stt_solve1
//(
//	double*x			// inout:	[-]:	解きたい関数の引数．
//) {
//	// ソルバの初期化．
//	_TRNSP_HANDLE_t handle;
//	dtrnlsp_init(&handle, &stt.set[1].n, &stt.set[1].m, x, stt.set[1].eps, &stt.set[1].iter1, &stt.set[1].iter2, &stt.set[1].rs);
//
//	// 収束計算のループ．関数値およびヤコビアン配列の動的確保．
//	double *fvec = new double[stt.set[1].m];
//	double *fjac = new double[stt.set[1].m * stt.set[1].n];
//	int RCI_Request = 0;
//
//	while (true) {
//		dtrnlsp_solve(&handle, fvec, fjac, &RCI_Request);
//
//		if (RCI_Request == -1) {
//			printf("最大繰り返し数 %d を超過しました．収束に失敗したためプログラムを終了します．\n", stt.set[1].iter1);
//			break;
//		}
//		if (RCI_Request == -2) {
//			printf("修正量が %e [m/s] を下回りました．収束に失敗したためプログラムを終了します．\n", stt.set[1].eps[0] * Rigid::l / Rigid::t);
//			break;
//		}
//		if (RCI_Request == -3) {
//			printf("部材の力の釣り合いの平均が基準値 %f ％ を下回りました．計算は収束しました．\n", 100 * sqrt(stt.set[1].eps[1]) / (stt.set[1].m));
//			break;
//		}
//		if (RCI_Request == -4) {
//			printf("ヤコビアンの行列式が %e [N] を下回りました．収束に失敗したためプログラムを終了します．\n", sqrt(stt.set[1].eps[2]));
//			break;
//		}
//		if (RCI_Request == -5) {
//			printf("探索域が %e [m/s] を下回りました．収束に失敗したためプログラムを終了します．\n", sqrt(stt.set[1].eps[3])* Rigid::l / Rigid::t);
//			break;
//		}
//		if (RCI_Request == -6) {
//			printf("1ステップの変化量が %e [N] を下回りました．収束に失敗したためプログラムを終了します．\n", sqrt(stt.set[1].eps[4]));
//			break;
//		}
//		if (RCI_Request == 1)
//			Stt_Eq1(&stt.set[1].m, &stt.set[1].n, x, fvec);
//
//		if (RCI_Request == 2)
//			djacobi(Stt_Eq1, &stt.set[1].n, &stt.set[1].m, fjac, x, &stt.set[1].jac_eps);
//	}
//
//	// ソルバのシャットダウン．
//	dtrnlsp_delete(&handle);
//
//	// 配列の開放．
//	delete[] fvec, fjac;
//	MKL_Free_Buffers();
//
//	return RCI_Request;
//}
// MKL "dtrnlsp_solve" を使った静解析ソルバ．使い方はマニュアル参照．（変数のアサインをマニュアルと統一させています．）
//int BS_Calculator::Stt_solve2
//(
//	double*x			// inout:	[-]:	解きたい関数の引数．
//) {
//	// ソルバの初期化．
//	_TRNSP_HANDLE_t handle;
//	dtrnlsp_init(&handle, &stt.set[2].n, &stt.set[2].m, x, stt.set[2].eps, &stt.set[2].iter1, &stt.set[2].iter2, &stt.set[2].rs);
//
//	// 収束計算のループ．関数値およびヤコビアン配列の動的確保．
//	double *fvec = new double[stt.set[2].m];
//	double *fjac = new double[stt.set[2].m * stt.set[2].n];
//	int RCI_Request = 0;
//
//	while (true) {
//		dtrnlsp_solve(&handle, fvec, fjac, &RCI_Request);
//
//		if (RCI_Request == -1) {
//			printf("最大繰り返し数 %d を超過しました．収束に失敗したためプログラムを終了します．\n", stt.set[2].iter1);
//			break;
//		}
//		if (RCI_Request == -2) {
//			printf("修正量が %e [m/s] を下回りました．収束に失敗したためプログラムを終了します．\n", stt.set[2].eps[0] * Rigid::l / Rigid::t);
//			break;
//		}
//		if (RCI_Request == -3) {
//			printf("部材の力の釣り合いが基準値 %e [N] を下回りました．計算は収束しました．\n", sqrt(stt.set[2].eps[1]));
//			break;
//		}
//		if (RCI_Request == -4) {
//			printf("ヤコビアンの行列式が %e [N] を下回りました．収束に失敗したためプログラムを終了します．\n", sqrt(stt.set[2].eps[2]));
//			break;
//		}
//		if (RCI_Request == -5) {
//			printf("探索域が %e [m/s] を下回りました．収束に失敗したためプログラムを終了します．\n", sqrt(stt.set[2].eps[3])* Rigid::l / Rigid::t);
//			break;
//		}
//		if (RCI_Request == -6) {
//			printf("1ステップの変化量が %e [N] を下回りました．収束に失敗したためプログラムを終了します．\n", sqrt(stt.set[2].eps[4]));
//			break;
//		}
//		if (RCI_Request == 1)
//			Stt_Eq2(&stt.set[2].m, &stt.set[2].n, x, fvec);
//
//		if (RCI_Request == 2)
//			djacobi(Stt_Eq2, &stt.set[2].n, &stt.set[2].m, fjac, x, &stt.set[2].jac_eps);
//	}
//
//	// ソルバのシャットダウン．
//	dtrnlsp_delete(&handle);
//
//	// 配列の開放．
//	delete[] fvec, fjac;
//	MKL_Free_Buffers();
//
//	return RCI_Request;
//}