/******************************************************************************* 
!								"B4P_BearingCalculator.cpp"	
!													2018/11/27	[Core-T]	楠崎
!
!	4点接触玉軸受の計算をするオブジェクト．
!	ここで指す”calculate”とは接触計算や摩擦計算などの直接的な物理現象ではなく，
!	より高次の常微分方程式の解など数学的な演算と定義している．
!	このクラス全体のメイン関数に相当するStt_solve()とDyn_solve()では
!	下位クラスであるBearingObjectと外部荷重Fと拘束条件y=y0を入力した際の解を導出．
!	※本オブジェクトの役割は一般的な連立方程式f(x)=bや連立微分方程式df/dx=0などの解の導出と
!	　導出した解のBearingオブジェクトへの代入に限定．
!	　方程式の左辺f(x)やdf/dx，ヤコビアンなどの具体的な立式は下位のクラスで実施．
!	※方程式の拘束条件はこのクラスで定義
!*******************************************************************************/

#include "B4P_Calculator.h"

// このオブジェクト内でのグローバル変数的役割．
B4P_Bearing B4P_Calculator::b4p = B4P_Bearing();



//// 静解析を行う．
//void B4P_Calculator::Stt_solve(double tr, int lim, double dx, double dv, double dq, double dw) {
//	Stt_calc_stf(tr, lim, dx, dv, dq, dw);// STEP1:剛性計算
//	//Stt_calc_frc(tr, lim, dx, dv, dq, dw);// STEP2:摩擦計算
//	//Stt_calc_stf_frc(tr, lim, dx, dv, dq, dw);// STEP3:両方計算
//	return;
//}

//// 解を求めたい関数 F(X) を表現
//// int    *m, 		// in:	[-]:	出力配列サイズ．
//// int    *n, 		// in:	[-]:	入力配列サイズ．
//// double *x,		// in:	[-]:	関数の引数．
//// double *f		// out:	[-]:	関数の戻り値．
//void B4P_Calculator::func (MKL_INT *m, MKL_INT *n, double *x, double *f){
//	b4p.set_param_stf(x);
//	b4p.get_F_stf(f);
//	return;
//}

void B4P_Calculator::Stt_init(const B4P_FileIn::Static*SttSet) {
	// 静解析設定
	for (int i = 0; i < 3; i++) {
		// ソルバの設定（長さ6の配列）．各要素の意味はマニュアル参照．
		this->SttSet[i].eps[0] = SttSet[i].eps[0];
		this->SttSet[i].eps[1] = Numeric::Square(SttSet[i].eps[1]);
		this->SttSet[i].eps[2] = Numeric::Square(SttSet[i].eps[2]);
		this->SttSet[i].eps[3] = Numeric::Square(SttSet[i].eps[3]);
		this->SttSet[i].eps[4] = Numeric::Square(SttSet[i].eps[4]);
		this->SttSet[i].eps[5] = SttSet[i].eps[5];

		// 最大繰り返し数．
		this->SttSet[i].iter1 = SttSet[i].iter1;

		// 最大試行回数．
		this->SttSet[i].iter2 = SttSet[i].iter2;

		// 初期ステップ条件．
		this->SttSet[i].rs = SttSet[i].rs;

		// ヤコビアン計算の精度．
		this->SttSet[i].jac_eps = SttSet[i].jac_eps;

		// 配列サイズ．
		this->SttSet[i].n = SttSet[i].nX;
		this->SttSet[i].m = SttSet[i].nX;
	}
	return;
}

// 静解析(1) 剛性計算：F(X) の数値を求める．
void B4P_Calculator::Stt_Eq_stf(
	int    *m, 		// in:	[-]:	出力配列サイズ．
	int    *n, 		// in:	[-]:	入力配列サイズ．
	double *x,		// in:	[-]:	変数 X
	double *f		// out:	[-]:	関数 F(X) の計算値．
) {
	B4P_Calculator::b4p.set_Xstf(x);	// 変数 x の入力
	B4P_Calculator::b4p.get_Fstf(f);	// 関数 F(X) の計算
	return;
}

// 静解析 剛性計算
int B4P_Calculator::Stt_solve_stf(
	double*x			// inout:	[-]:	解きたい関数の引数．	
) {
	// ソルバの初期化．
	_TRNSP_HANDLE_t handle;
	MKL_INT res = dtrnlsp_init(&handle, &SttSet[0].n, &SttSet[0].m, x, SttSet[0].eps, &SttSet[0].iter1, &SttSet[0].iter2, &SttSet[0].rs);
	if(res != TR_SUCCESS)
	{
	cout << "error in dtrnlsp_initn";
	cin.get();
	return 0;
	}
	// 収束計算のループ．関数値およびヤコビアン配列の動的確保．
	double *fvec = new double[SttSet[0].m];
	double *fjac = new double[SttSet[0].m * SttSet[0].n];
	int RCI_Request = 0;

	while (true) {
		dtrnlsp_solve(&handle, fvec, fjac, &RCI_Request);

		if (RCI_Request == -1) {
			printf("最大繰り返し数 %d を超過しました．収束に失敗したためプログラムを終了します．\n", SttSet[0].iter1);
			break;
		}
		if (RCI_Request == -2) {
			printf("修正量が %e [μm] を下回りました．収束に失敗したためプログラムを終了します．\n", SttSet[0].eps[0]);
			break;
		}
		if (RCI_Request == -3) {
			printf("部材の力の釣り合いが基準値 %e [N] を下回りました．計算は収束しました．\n", sqrt(SttSet[0].eps[1]));
			break;
		}
		if (RCI_Request == -4) {
			printf("ヤコビアンの行列式が %e [N] を下回りました．収束に失敗したためプログラムを終了します．\n", sqrt(SttSet[0].eps[2]));
			break;
		}
		if (RCI_Request == -5) {
			printf("探索域が %e [μm] を下回りました．収束に失敗したためプログラムを終了します．\n", sqrt(SttSet[0].eps[3]));
			break;
		}
		if (RCI_Request == -6) {
			printf("1ステップの変化量が %e [N] を下回りました．収束に失敗したためプログラムを終了します．\n", sqrt(SttSet[0].eps[4]));
			break;
		}
		if (RCI_Request == 1)
			Stt_Eq_stf(&SttSet[0].m, &SttSet[0].n, x, fvec);

		if (RCI_Request == 2)
			djacobi(Stt_Eq_stf, &SttSet[0].n, &SttSet[0].m, fjac, x, &SttSet[0].jac_eps);
	}

	// ソルバのシャットダウン．
	dtrnlsp_delete(&handle);

	// 配列の開放．
	delete[] fvec, fjac;
	MKL_Free_Buffers();

	return RCI_Request;
}



// 動解析の計算条件を示す定数を設定
void B4P_Calculator::Dyn_init(const B4P_FileIn::Dynamic&DynSet) {

	this->dyn_h		= DynSet.h;		// 最小ステップサイズ
	this->dyn_hmin	= DynSet.hmin;	// 許容誤差
	this->dyn_ep	= DynSet.ep;	// しきい値
	this->dyn_tr	= DynSet.tr;	// 初期ステップサイズ
	this->dyn_n = DynSet.nX;		// 変数yの要素長
	for (int i = 0; i < 128; i++)
		this->ipar[i] = 0;
	this->ipar[1] = 1;
	this->dpar = new double[13 * this->dyn_n];

	this->t_step = DynSet.calctime / (double)DynSet.sampling;	// サンプリング時間ステップサイズ（≠最小時間ステップサイズ）
	this->calctime = DynSet.calctime;	// 計算区間（時間）
	return;
}


//// Intel ODE Solver を用いて動解析を行う．（内部の記載はほとんどコピペ．dodesolは取得後変更無し．）
void B4P_Calculator::Dyn_solve(double*y,double*t) {
	// 以下の変数のアサインはMKLマニュアルに準ずる．詳細はマニュアルを参照のこと．
	int    ierr, kd[2];
	double t0 = t[0] / Rigid::t;
	double t1 = t[1] / Rigid::t;
	dodesol(this->ipar,&this->dyn_n, &t0,&t1,y,&B4P_Calculator::Dyn_Eq,&B4P_Calculator::Dyn_void,
		&this->dyn_h,&this->dyn_hmin,&this->dyn_ep,&this->dyn_tr,this->dpar,kd,&ierr);

	t[0] = t0 * Rigid::t;
	t[1] = t1 * Rigid::t;
	return;
}


// ベクトル表現による運動方程式の定義．12*2+9*n次元で構成されている．
// intel_odeライブラリの書式に準じているため，詳細はそちらのマニュアル参照．
void B4P_Calculator::Dyn_Eq(int*n,double*t0,double*y,double*dydt) {

	// まず入力値から各部材の変位を更新．
	b4p.set_y(y);

	// 更新した変位から荷重を計算．
	b4p.get_dydt(dydt);
}


// ダミー関数．引数に必要だが使わないため適当に定義した．
void B4P_Calculator::Dyn_void(void) {
}

B4P_Calculator::~B4P_Calculator() {
	delete[] dpar;
}









































//// このメソッドは一般的な方程式f(x)<0を解く単純なメソッドに置き換える
//void B4P_Calculator::Stt_calc_approx(double *y_stf, int lim){
//	//double*F_stf = new double[b4p.nX_stf];
//	double appro_dx = 1.0e-6;
//
//	// ヤコビアン試算用（後で消す）
//	double*Jacobian = new double[b4p.nX_stf * b4p.nX_stf];
//	MatrixXd Jacobian_(b4p.nX_stf, b4p.nX_stf);
//	double tr = 1e-3;	//誤差閾値
//	double dx  = 1e-9;	//変位微小変化量
//	double dv = 1e-9;	//速度微小変化量
//	double dq = 1e-9;	//クォータニオン微小変化量
//	double dw = 1e-9;	//角速度微小変化量
//	double*F_stf = new double[b4p.nX_stf]; 
//	VectorXd F_stf_(b4p.nX_stf);
//	double *F_ball = new double[3];
//	// 近似解1-1.遠心力・接触力の釣り合いの式から玉径方向座標を概算
//
//	// ball[0]と外輪が接触し，(遠心力)<(接触力)となるまでball[0]を移動
//	//（回転速度が0の場合や，玉中心が外輪軌道よりも外側にある場合は想定せず）
//	// 一定以上のループを繰り返したら計算を中止するようにしたほうがいいかもしれない．
//	while(true){
//		b4p.set_param_stf(y_stf);
//		//ball_fn = b4p.get_fn_ball();
//		b4p.get_F_stf_ball(F_ball,0);
//		//cout << y_stf[6] <<":"<< ball_fn <<endl;
//		if(F_ball[5] < 0)
//			break;
//		y_stf[6]+= appro_dx;
//	}
//
//
//
//
//	// 近似解1-2.他の玉に対してもball[0]と同じ座標を代入
//	for(int i = 1; i < b4p.Z; i++)
//		y_stf[2 * i + 6] = y_stf[6];
//
//
//	// 近似解1-3.内輪の力の釣り合いの式から内輪yz座標を求める．
//	// （ただし，Y方向またはZ方向の荷重がないときはスキップ）
//	//	F(荷重)-N(ボールとの接触力)=0
//	// 一定以上のループを繰り返したら計算を中止するようにしたほうがいいかもしれない．
//	double fn_ir;
//	// 内輪の外部荷重方向の単位ベクトルを取得
//	Vector3d f_dir = b4p.get_Fyzload_dir();
//	if(f_dir[1] !=0 || f_dir[2] != 0){
//		while(true){
//			b4p.set_param_stf(y_stf);
//			// 内輪にかかる接触力・外部荷重の合力を計算
//			b4p.get_F_stf(F_stf);
//			fn_ir = Vector3d(0,F_stf[1],F_stf[2]).dot(f_dir);
//			// 外部荷重より接触力が大きくなったら計算終了
//			if(fn_ir < 0)
//				break;
//			// 内輪を外部荷重の向きに移動
//			y_stf[1]+= appro_dx*f_dir[1];
//			y_stf[2]+= appro_dx*f_dir[2];
//		}
//	}
//
//	// 近似解2.玉・内輪x座標を求める．
//	// ただし，X方向の荷重がないときはスキップ
//	double fx_ball = 0;
//	int fx_dir = b4p.get_Fxload_dir();
//	if(fx_dir != 0){
//		// 近似解2-1.玉を移動させる
//		while(true){
//			b4p.set_param_stf(y_stf);
//			// 玉にかかる接触力・外部荷重の合力を計算
//			b4p.get_F_stf(F_stf);
//			for(int i = 0; i < b4p.Z; i++){
//				fx_ball += F_stf[2*i+5];
//			}
//			//cout << y_stf[5] << "," << y_stf[7] << ":"<< fx_ball <<endl;
//			// 外部荷重より外輪との接触力が大きくなったら計算終了
//			if(fx_ball < 0) 
//				break;
//			// すべての玉を外部荷重の向きに移動
//			for(int i = 0; i < b4p.Z; i++)
//				y_stf[2 * i + 5] += appro_dx * fx_dir;
//		}
//		// 近似解2-2.内輪を移動させる
//		while(true){
//			b4p.set_param_stf(y_stf);
//			// 内輪にかかる接触力・外部荷重の合力を計算
//			b4p.get_F_stf(F_stf);
//			cout << F_stf[0] << F_stf[1] << F_stf[2] << endl;
//			cout << y_stf[0] << y_stf[1] << y_stf[2] << endl;
//			// 外部荷重より玉との接触力が大きくなったら計算終了
//			if(F_stf[0] < 0) 
//				break;
//			// 内輪を外部荷重の向きに移動
//			y_stf[0]+= appro_dx * fx_dir;
//		}
//	}
//	//近似解3. 内輪の姿勢を求める．
//	// 
//	Vector3d myz_dir = b4p.get_nyzload_dir();
//	if(myz_dir[1] !=0 || myz_dir[2] != 0){
//		Vector3d myz_ir;
//		double ay =0;
//		double az =0;
//		for(int i = 0; i < lim; i++){
//			b4p.set_param_stf(y_stf);
//			// 内輪にかかる接触モーメント・外部荷重の合力を計算
//			//myz_ir =  Vector2d(F_stf[3], F_stf[4]); //b4p.get_myz_ir();
//			//b4p.get_F_stf(F_stf);
//			double M_ = Vector3d(0, F_stf[3], F_stf[4]).dot(myz_dir);
//			//cout << ay <<", " << az <<", "<< y_stf[3] << ", " << y_stf[4] << ", "<< myz_ir[1]  << ", " <<  myz_ir[2] << ", " << M_ <<endl;
//			//cout << y_stf[3] << ", " << y_stf[4] << ", "<< myz_ir[0]  << ", " <<  myz_ir[1] <<endl;
//			// 外部荷重より玉との接触力が大きくなったら計算終了
//			if(myz_ir[1] < 0 && myz_ir[2] < 0)
//				break;
//			// 内輪の姿勢を外部荷重の向きに移動
//			ay += appro_dx * pow(2, i) * myz_dir[2];
//			az += appro_dx * pow(2, i) * -myz_dir[1];
//			b4p.set_q_ir(Vector3d (1, ay, az));
//			b4p.get_param_stf(y_stf);
//		}
//
//	}
//	return;
//}

//// 静解析(STEP1:剛性計算)を行う．
//void B4P_Calculator::Stt_calc_stf(double tr, int lim, double dx, double dv, double dq, double dw) {
//
//
//
//	// まずは初期値として現在の姿勢の配列を取得
//	double*y_stf = new double[b4p.nX_stf];
//	b4p.get_param_stf(y_stf);
//
//
//
//
//
//	// 変数の確保．
//	VectorXd dy_stf_(b4p.nX_stf);
//	double*F_stf = new double[b4p.nX_stf]; 
//	VectorXd F_stf_(b4p.nX_stf);
//	double*Jacobian = new double[b4p.nX_stf * b4p.nX_stf];
//	MatrixXd Jacobian_(b4p.nX_stf, b4p.nX_stf);
//	VectorXd y_stf_(b4p.nX_stf);		// 後で絶対に消す
//	VectorXd dy_stf_2(b4p.nX_stf);		// 1回前の修正量
//	dy_stf_2.setZero();
//
//	// y_stfが解の近傍になるように近似値を計算
//	Stt_calc_approx(y_stf, lim);
//	
//	double *x = new double [b4p.nX_stf];
//	double *LW = new double [b4p.nX_stf];
//	double *UP = new double [b4p.nX_stf];
//	double eps[6];
//	for(int i = 0; i < b4p.nX_stf; i++){
//		x[i] = y_stf_[i];
//		LW[i] = -1000;
//		UP[i] = 1000;
//	}
//	eps[0] = 0.0001;
//	eps[1] = 0.0001;
//	eps[2] = 0.0001;
//	eps[3] = 0.0001;
//	eps[4] = 0.0001;
//	eps[5] = 0.0001;
//	MKL_INT iter1 = 1000, iter2 = 100;	// 最大繰り返し数
//	double rs = 0.0;					// 初期ステップサイズ
//	double *fvec = NULL;				// F(X)の数値
//	double *fjac = NULL;				// ヤコビアン
//	MKL_INT RCI_Request;
//	int m = b4p.nX_stf, n = b4p.nX_stf;
//	// ソルバーの初期化
//	_TRNSPBC_HANDLE_t handle;			// ソルバー内部の情報が格納された構造体
//	MKL_INT init_res = dtrnlspbc_init (&handle, &b4p.nX_stf, &b4p.nX_stf, x, LW, UP, eps, &iter1, &iter2, &rs);
//	if (init_res != TR_SUCCESS){
//	    printf ("| error in dtrnlspbc_init\n");
//	    MKL_Free_Buffers ();
//	    return;
//	}
//
//	int successful = 0;
//	while (successful == 0)
//	{
//		MKL_INT sol_res = dtrnlspbc_solve (&handle, fvec, fjac, &RCI_Request);
//		// 途中経過確認用
//		for(int i = 0; i < n; i++){
//			cout  << "  x" << i << "= " << x[i];
//		}
//		cout << endl;
//
//		if (sol_res != TR_SUCCESS){
//	        printf ("| error in dtrnlspbc_solve\n");  
//	        MKL_Free_Buffers ();
//	        /* and exit */
//	        return;
//	    }
//	    /* according with rci_request value we do next step */
//	    if (RCI_Request == -1 ||
//	        RCI_Request == -2 ||
//	        RCI_Request == -3 ||
//	        RCI_Request == -4 || RCI_Request == -5 || RCI_Request == -6){
//	        successful = 1;
//		}
//		if (RCI_Request == 1){
//	        func (&m, &n, x, fvec);
//	    }
//	    if (RCI_Request == 2){
//			MKL_INT jac_res = djacobi (func, &n, &m, fjac, x, eps);
//	        if (jac_res != TR_SUCCESS)
//	        {
//	            printf ("| error in djacobi\n");
//	            MKL_Free_Buffers ();
//	            return;
//	        }
//	    }
//
//
//	}
//
//
//	// 後で絶対に消す．
//	ofstream file;
//	string name = "Stt.csv";			// [in]	ファイル名
//	file.open(name.c_str(), ios_base::trunc);
//	string head = "\n";
//	file << head;
//
//	// 収束するまでy_stfを調整する．
//	for (int loop = 0; loop < lim; loop++) {
//		// 現在のy_stfに対してf(y_stf)の数値を導出
//		b4p.set_param_stf(y_stf);
//		b4p.get_F_stf(F_stf);
//		for (int i = 0; i < b4p.nX_stf; i++)
//			F_stf_(i) = F_stf[i];
//
//
//		cout << "loop = " << loop << "	";
//		cout << F_stf_.norm() << endl;
//		// 収束（荷重配列のノルムがしきい値以下）していたら，ループを抜ける．
//		if (F_stf_.norm() < tr)
//			break;
//		// ヤコビアン行列J(f(y))を導出
//		b4p.get_Jacobian_stf(y_stf, dx, dq, Jacobian);
//
//		// 逆行列計算のため，double配列からEigen::MatrixXdにコピー
//		for (int i = 0; i < b4p.nX_stf; i++)
//			for (int j = 0; j < b4p.nX_stf; j++)
//				Jacobian_(i, j) = Jacobian[b4p.nX_stf * i + j];
//
//		// 数値確認用．後で絶対消す．
//		VectorXd vector(b4p.nX_stf*2);	// [in]	書き込む数値配列
//		int prec = 10;				// [in]	数値の有効数字
//		for (int i = 0; i < b4p.nX_stf; i++)
//			y_stf_(i) = y_stf[i];
//		vector << y_stf_, F_stf_;
//		IOFormat CSVFormat(prec, 0, ", ", ", ", "", "", "", "\n");
//		ofstream file;
//		file.open(name.c_str(), ios_base::app);
//		file << vector.format(CSVFormat);
//
//		// ヤコビアンの逆行列が導出できない場合，エラーを出力し計算を終了
//		if (Jacobian_.determinant() == 0){
//			cout << Jacobian_.determinant()  << endl;
//			cout << "static analysis calcation error!" << endl;
//
//			cout << "____________________" << endl;
//			break;
//		}
//
//		dy_stf_ = Jacobian_.colPivHouseholderQr().solve(F_stf_);
//		dy_stf_(4) = 0;
//		dy_stf_(3) = 0;
//		dy_stf_(2) = 0;
//		dy_stf_(1) = 0;
//		for (int i = 0; i < b4p.nX_stf; i++)
//			y_stf[i] -= dy_stf_(i)*0.3; //- dy_stf_2(i)*0.15;
//
//		for (int i = 0; i < b4p.nX_stf; i++)
//			dy_stf_2(i) = dy_stf_(i);
//	}
//}
//
