/*******************************************************************************
"BS_main.cpp"
2018/12/05	[Core-T]	楠崎

4点接触玉軸受の静・動解析を行うプログラム．

!*******************************************************************************/

/*******************************************************************************
"ヘッダーファイル" の使い方
2019/04/15	 [Core-T]	楠崎

1. 右クリック"プロパティ" → "追加のインクルードディレクトリ"
2. "\（現プロジェクト）\include"

!*******************************************************************************/

/*******************************************************************************
"Eigen3.3.7" の使い方
2019/04/15	 [Core-T]	楠崎

1. VersionManagerから"Eigen"を任意の階層に"取得"
2. 右クリック"プロパティ" → "追加のインクルードディレクトリ"
3. "\（取得した階層）\Engen"

※1. Eigenオブジェクトは値渡しができません．"const Vector3d&x"のように参照渡しで記述して下さい．
※2. Eigenをメンバにもつオブジェクトには"EIGEN_MAKE_ALIGNED_OPERATOR_NEW"のマクロを定義して下さい．

!*******************************************************************************/

/*******************************************************************************
"intel_ode" の使い方
2019/04/15	 [Core-T]	楠崎

1. VersionManagerから"intel_ode"を任意の階層に"取得"
2. 以下の手順に従う．
（2019/04/15　Intel公式ホームページから引用）
インテル® マス・カーネル・ライブラリー (インテル® MKL) をリンクするために Microsoft* Visual C/C++* 開発システムを設定する手順は、インテル® Parallel Studio XE Composer Edition の Microsoft* Visual Studio* 統合コンポーネントをインストールしたかどうかに依存します。
統合コンポーネントをインストールした場合は、「Microsoft* Visual C/C++* プロジェクトとインテル® MKL の自動リンク」を参照してください。
統合コンポーネントをインストールしていない場合、またはインテル® MKL をリンクするために別途設定が必要な場合は、次の手順を実行して、Microsoft* Visual C++* 開発システムの設定を行います。Visual C++* のバージョンによって以下で説明している一部のメニュー項目は異なりますが、基本的な設定手順はすべてのバージョンで同じです。
[ソリューション エクスプローラー] でプロジェクトを右クリックして [プロパティ] をクリックします。
[構成プロパティ] > [VC++ ディレクトリ] を選択します。
[インクルード ディレクトリ] を選択します。 インテル® MKL インクルード・ファイルのディレクトリー (<mkl ディレクトリー>\include) を追加します。
[ライブラリ ディレクトリ] を選択します。 インテル® MKL ライブラリーと OpenMP* ライブラリーのアーキテクチャー固有のディレクトリー (例えば、<mkl ディレクトリー>\lib\ia32 および <親製品のディレクトリー>\compiler\lib\ia32) を追加します。
[実行可能ファイル ディレクトリ] を選択します。 ダイナミック・リンク・ライブラリーのアーキテクチャー固有のディレクトリーを追加します。
OpenMP* サポートの場合、例えば <親製品のディレクトリー>\redist\ia32\compiler と入力します。
インテル® MKL の場合 (ダイナミックにリンクする場合のみ)、例えば <親製品のディレクトリー>\redist\ia32\mkl と入力します。
[構成プロパティ] > [Custom Build Step (カスタム・ビルド・ステップ)] > [Additional Dependencies (追加の依存ファイル)] を選択します。 必要なライブラリー (例えば、mkl_intel_c.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib) を追加します。

!*******************************************************************************/

/*******************************************************************************
おすすめ設定
2019/04/15	 [Core-T]	楠崎

1. ソリューションディレクトリに \out フォルダを作成．
2. 作業ディレクトリを "$(ProjectDir)out" に設定．
3. 出力ディレクトリを "$(ProjectDir)out\$(Configuration)\" に設定
4. 中間ディレクトリを "$(ProjectDir)out\$(Configuration)\$(ProjectName)\" に設定
5. 計算に必要な入力ファイルは全て \out 直下に入れましょう．

C:\Program Files (x86)\Microsoft Visual Studio 11.0\Common7\IDE に テキストファイルから usertype.dat を作成し，その中に
Vector2d, Matrix2d, Array2d, Vector3d, Matrix3d, Array3d, Vector4d, Matrix4d, Array4d, VectorXd, MatrixXd, ArrayXd, Quaterniond, string
をカンマなし・改行で記入．
!*******************************************************************************/

/*******************************************************************************
Visual Studio2017 にする際の注意点１
2019/12/04	 [Core-T]	楠崎

Windows8.1 SDK が足りないと文句を言われます．以下の方法で追加してください．（参考：http://tooljp.com/qa/MSB8036-build-compile-error-88E1.html）

Windows8.1 SDK を追加でインストールします。

(1)Visual Studio 2017のメニューから [ファイル] - [新規作成] - [プロジェクト] を選択します。
(2)[Visual Studio インストーラを開く]を選択します。
(3)[個別コンポーネント]タブを選択します。
(4)[Windows 8.1SDK]のチェックをオンにします。
!*******************************************************************************/

/*******************************************************************************
その他注意事項
2019/04/15	 [Core-T]	楠崎

C言語ライブラリの使用でLNK2019エラーが出る場合はextern"C"を参考にして下さい．
自分は全部それでした．

ループが多用される箇所に "VectorXd" を使わないこと！めっちゃ遅くなります！

!*******************************************************************************/

#define _MAX_STRING_ 256
#define _DYN_LAST_   1

#include <iostream>			// cmd画面出力のため．
#include <string.h>			// ファイル名操作のため．
#include <string>			// ファイル名操作のため．
#include <conio.h>			// Enterキー入力でプログラム実行を止めるため．
#include <chrono>			// 実行時間計測のため．
#include <Eigen\Dense>		// Tempファイル出力のため．

using Eigen::VectorXd;		// Tempファイル出力のため．

#include "BS_FileIn.h"
#include "BS_FileOut.h"
#include "BS_Calculator.h"
#include "BS_BallScrew.h"
#include "Rigid.h"
#include "LicenceSimple.h"

int sub_Dynamic(BS_Calculator & calc, BS_FileOut&FO, int i, double*t, double*y, string&Temp, bool stopcalc, VectorXd&y_, double dTerr, int stp);

int main(int argc, char *argv[]) {

	// タイトル表示．
	printf("       ---------------------------------------------------\n");
	printf("\n");
	printf("            DYNAMIC BEARING ANALYSIS IN NSK (D-BRAIN)\n");
	printf("\n");
	printf("              D-BS    : %s (Build Date)\n", __DATE__);
	printf("\n");
	printf("            Dynamic motion analysis of Ball Screw \n");
	printf("\n");
	printf("       ---------------------------------------------------\n");

	if (!LicenceSimple::hasTime(2023, 5, 1))
		return -1;

	char fpath_dbsin_csv[_MAX_STRING_];
	if (argc < 2) {
		std::cout << "インプットファイルの名前を入力してください" << std::endl;
		cin >> fpath_dbsin_csv;
	}
	else
		sprintf_s(fpath_dbsin_csv, _MAX_STRING_, argv[1]);

	BS_FileIn FI;

	//CSVファイルから数値を読み込み，FIの各メンバ変数に書き込む
	try {
		FI.read_input_all(fpath_dbsin_csv);
	}
	catch (exception&e) {
		std::cout << e.what() << std::endl;
		std::cout << "inputファイル読み取り時にエラーが発生したため強制終了します．" << std::endl;
		return -1;
	}

	// 世界全体の単位系の統一．
	Rigid::l = FI.rigid.l;
	Rigid::t = FI.rigid.t;
	Rigid::g = Vector3d(FI.rigid.g);

	BS_Calculator calc;

	// ボールねじオブジェクトの初期化．
	calc.BS.allocate(FI);
	calc.BS.init(FI, FI.stt.v0, FI.stt.w0, FI.stt.wn);

	// 出力形式について初期化．
	BS_FileOut FO;
	FO.init(FI);
	BS_Calculator::init_stt(FI.stt, FI.ballnum);
	BS_Calculator::init_dyn(FI.dyn, FI.ballnum);

	// シャフト初期位置決定．
	string Temp = FI.output.temp + ".csv";
	double t_str = 0; // 初期時刻[s]
	std::cout << std::endl << "$$PositionSet の設定が" << int(FI.initial.preset) << "であったため";
	switch (FI.initial.preset) {
	case BS_FileIn::Initial::Preset::ReadPos:
		std::cout << "シャフト初期位置を入力値へ移動します．" << std::endl;
		calc.BS.lock_y0(FI.initial.x0, FI.initial.ax0, FI.stt.v0, FI.stt.w0);
		break;
	case BS_FileIn::Initial::Preset::ReadTemp:
		std::cout << "tempファイルから玉・シャフトの位置・速度と現在時刻を取得します．" << std::endl;
		double  *_y = new double[BS_Calculator::dyn.set[_DYN_LAST_].nX];
		BS_FileOut::read_Params(Temp, BS_Calculator::dyn.set[_DYN_LAST_].nX, t_str, _y);
		calc.BS.set_dyn_y1(_y);
		break;
	}

	// 静解析．
	std::cout << std::endl << "【D-BS 静解析】" << std::endl;

	if (!FI.runStatic)
		std::cout << "は $$SttMode = 0 であったため行いません．" << std::endl;

	else {
		// step0（玉一様移動簡易計算）
		std::cout << std::endl << "【Approximate0：玉一様移動簡易計算】" << std::endl;
		if (!FI.stt.run[0])
			std::cout << "は入力が0であったため行いません．" << std::endl;

		else {
			std::cout << "おおよその位置へシャフトを移動します．" << std::endl;

			calc.BS.preset_y0(1e-9, 1e-12, 1e-9);
		}
		// step1（剛性計算）．MKL "dtrnlsp_solve" を使って解く．
		std::cout << std::endl << "【Statics0：摩擦なし剛性計算】" << std::endl;
		if (!FI.stt.run[1])
			std::cout << "は入力が0であったため行いません．" << std::endl;

		else {
			double *x0 = new double[BS_Calculator::stt.set[0].n];
			calc.BS.get_y0(x0);
			int RCI_Request = BS_Calculator::Stt_solve(0, x0);
			delete[] x0;
			if (RCI_Request != -3)
				if (!FI.stt.run[2])
					return -1;
				else
					std::cout << "収束はしませんでしたが次のSTEP1で釣り合う可能性もあるため，計算を続行します．" << std::endl;
		}
		// step2（摩擦計算）．MKL "dtrnlsp_solve" を使って解く．
		std::cout << std::endl << "【Dynamics0：簡易摩擦計算】" << std::endl;
		if (!FI.stt.run[2])
			std::cout << "は入力が0であったため行いません．" << std::endl;

		else {
			const int i = 0;
			double  *y = new double[BS_Calculator::dyn.set[i].nX];
			VectorXd y_temp = VectorXd(1);
			double t0[2];
			t0[0] = 0;	t0[1] = BS_Calculator::dyn.set[i].t_step;
			calc.BS.init_dyn0();
			calc.BS.get_dyn_y0(y);
			sub_Dynamic(calc, FO, i, t0, y, Temp, FI.dyn.set[i].stopcalc, y_temp, FI.dyn.set[i].dTerr, FI.dyn.set[i].stp);
			calc.BS.deinit_dyn0(FI.stt.v0, FI.stt.w0);
		}
		// step2.5（純転がり速度計算）．恒等式から各玉の主荷重2点での純転がり速度を求める．
		std::cout << std::endl << "【Approximate1：純転がり速度計算】" << std::endl;
		if (!FI.stt.run[3])
			std::cout << "は入力が0であったため行いません．" << std::endl;

		else {
			calc.BS.pure_Rolling();
			std::cout << "玉は全て純転がりに設定されました．" << std::endl;
		}
		// step3（定常状態）．恒等式から各玉の主荷重2点での純転がり速度を求める．
		std::cout << std::endl << "【Statics1：各玉定常状態計算】" << std::endl;
		if (!FI.stt.run[4])
			std::cout << "は入力が0であったため行いません．" << std::endl;

		else {
			if (!FI.stt.run[3])
				std::cout << "注意！純転がり計算をしないと収束しない危険性が高まります！設定を見直してください！" << std::endl;
			double *x1 = new double[BS_Calculator::stt.set[1].n];

			bool has_next = true; int ib = 0;
			while (has_next) {
				std::cout << "玉番号:\t" << ib << ",\t";
				has_next = calc.BS.get_y1(ib, x1);
				BS_Calculator::stt.i1 = ib;
				int RCI_Request = BS_Calculator::Stt_solve(1, x1);
				//if (RCI_Request != -3)
				//	return -1;
				ib++;
			}
			delete[] x1;
		}
	}
	// 動解析．
	double  *y = new double[BS_Calculator::dyn.set[_DYN_LAST_].nX];
	VectorXd y_ = VectorXd(BS_Calculator::dyn.set[_DYN_LAST_].nX + 1);
	double t1[2];
	t1[0] = t_str;
	t1[1] = BS_Calculator::dyn.set[_DYN_LAST_].t_step + t_str;
	y_[0] = t1[0];
	calc.BS.get_dyn_y1(y);

	// 前回計算結果を消去し，初期状態での座標・速度をtempに出力
	if (FI.output.deletelastout) {
		FileOut::write_header(Temp, "");
		for (int i = 0; i < BS_Calculator::dyn.set[_DYN_LAST_].nX; i++)
			y_[i + 1] = y[i];
		FileOut::write_vector(Temp, y_, 15);
	}

	std::cout << std::endl << "【D-BS 動解析】" << std::endl;
	sub_Dynamic(calc, FO, _DYN_LAST_, t1, y, Temp, FI.dyn.set[_DYN_LAST_].stopcalc, y_, FI.dyn.set[_DYN_LAST_].dTerr, FI.dyn.set[_DYN_LAST_].stp);

	std::cout << "ポスト処理に移行します．" << std::endl;

	// とりあえずベタ打ちでポスト処理．
	FO.write_AllHeader();
	ifstream ifs(Temp);
	string line;
	double*Params = new double[BS_Calculator::dyn.set[_DYN_LAST_].nX];
	while (getline(ifs, line)) {
		vector<string> strvec = FileIn::split(line, ',');
		for (int i = 0; i < BS_Calculator::dyn.set[_DYN_LAST_].nX; i++)
			Params[i] = stod(strvec.at(i + 1));
		calc.BS.set_dyn_y1(Params);

		calc.BS.get_dyn_dydt1(y, 0, 0);	// 0/0 は v/w の加減速パラメタ．ポスト処理に加減速は含まれないが，一応0で定常状態を出力させる．
		calc.BS.save(FO);
		FO.write_AllParams(stod(strvec.at(0)));
	}
	std::cout << std::endl << "全ての計算が完了しました．" << std::endl;

	// 動的確保した配列の解放．
	delete[] y, Params;

	// 正常終了．
	return 0;
}

// 動解析．
int sub_Dynamic(BS_Calculator & calc, BS_FileOut&FO, int i, double*t, double*y, string&Temp, bool stopcalc, VectorXd&y_, double dTerr, int stp) {

	std::cout << "\n指定ループ回数を達成するか，Enterキーを押せば計算を終了します．" << std::endl << std::endl;

	// 時間計測開始．
	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

	// 1回目の計算で打ち切りしないようにトルクの初期値は十分大きい数値に設定
	double Ti0 = 1e10, Ti1 = 0;
	int cnt = 0;

	while (t[0] < BS_Calculator::dyn.set[i].t_end) {

		bool finished = false;

		calc.Dyn_solve(y, t, i);

		// 次の時間ステップの設定．（t[0]はDyn_solve内で自動的にt[1]に変更される．）
		t[1] += BS_Calculator::dyn.set[i].t_step;
		std::cout << "\r" << "Now " << 1e3 * t[0] << " [ms].\t"
			<< int(100.0 * t[0] / BS_Calculator::dyn.set[i].t_end) << "% Complete." << std::string(20, ' ');

		// 1ステップ終了ごとにテキストに書き込む．
		if (i == _DYN_LAST_) {
			y_[0] = t[0];
			calc.BS.get_dyn_y1(y);

			for (int j = 0; j < BS_Calculator::dyn.set[i].nX; j++)
				y_[j + 1] = y[j];
			FileOut::write_vector(Temp, y_, 15);
		}
		// enterキーを押したらループ中断
		if (_kbhit() != 0 && _getch() == '\r')
			finished = true;

		// 計算自動終了がある場合，1ステップごとに収束判定
		if (stopcalc) {

			// 収束判定を行い，指定された回数連続で基準値を下回っていた場合はループ中断
			calc.BS.save(FO);
			Ti1 = FO.ST_CY[0].Ts[0];
			double dT = Ti1 - Ti0;
			cout << "dT = " << Unit::Nm2Nmm(dT) << " [Nmm]";

			cnt = (abs(dT) < dTerr) ? cnt + 1 : 0;

			std::cout << "\tcount = " << cnt;
			if (cnt >= stp) {
				std::cout << std::endl << "打ち切り基準値を" << stp
					<< "回連続で下回ったので計算を打ち切ります" << std::endl;
				finished = true;
			}
			Ti0 = Ti1;
		}
		if (finished)
			break;

	}
	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	double elapsed = double(std::chrono::duration_cast<std::chrono::seconds>(end - start).count());
	std::cout << std::endl << std::endl << "計算に " << elapsed << " [s] かかりました．" << std::endl << std::endl;

	return 0;
}


//const int i = 0;
//double  *y = new double[BS_Calculator::dyn.set[i].nX];
//VectorXd y_temp = VectorXd(1);
//double t0[2];
//t0[0] = 0;	t0[1] = BS_Calculator::dyn.set[i].t_step;
//calc.BS.init_dyn0();
//calc.BS.get_dyn_y0(y);
//sub_Dynamic(calc, FO, i, t0, y, Temp, FI.dyn.set[i].stopcalc, y_temp, FI.dyn.set[i].dTerr, FI.dyn.set[i].stp);
//calc.BS.deinit_dyn0();


//// step2（力の釣り合い計算）．MKL "dtrnlsp_solve" を使って解く．
//std::cout << std::endl << "【step2：力の釣り合い計算】" << std::endl;
//double *x2 = new double[calc.stt.set[2].n];
//for (int ib = 0; ib < n; ib++) {
//	std::cout << "Ball No." << ib << ":\t";
//	calc.BS.set_i1(ib);
//	calc.BS.get_y2(x2);
//	BS_Calculator::Stt_solve2(x2);
//}
//delete[] x2;
	//// 与圧入力の時だけ与圧変位の計算を行う．
	//if (FI.preload.mode == 1) {
	//	std::cout << std::endl << "【step p：与圧計算】" << std::endl;
	//	double *xp = new double[calc.sttsetp.n];
	//	calc.BS.get_yp(xp);
	//	BS_Calculator::Stt_solvep(xp);
	//	cout << "-x側ナットの変位は " << Unit::m2mm(xp[0]) << "[mm] です．" << endl
	//		<< "+x側ナットの変位は " << Unit::m2mm(xp[1]) << "[mm] です．" << endl;
	//	delete[] xp;
	//}


			//if (calc.BS.shaft_DOF() < 5)
			//	std::cout << "完全非拘束でない場合の簡易計算は推奨されません．ねじ位置は入力拘束値と異なる変位に固定されます．" << std::endl;

			//double *x1 = new double[BS_Calculator::stt.set[1].n];
			//calc.BS.get_y0(x1);
			//int RCI_Request = BS_Calculator::Stt_solve(1, x1);
			//delete[] x1;
			//if (RCI_Request != -3)
			//	return -1;


			//if (abs(dT) < dTerr) {
			//	cnt++;
			//}
			//else {
			//	cnt = 0;
			//}


