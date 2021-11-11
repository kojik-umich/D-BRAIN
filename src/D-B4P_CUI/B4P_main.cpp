/*******************************************************************************
"B4P_main.cpp"
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

#include <iostream>
#include <string.h>
#include <string>
#include <Eigen\Dense>
#include <conio.h>
#include <chrono>

#include "B4P_FileIn.h"
#include "B4P_FileOut.h"
#include "B4P_Calculator.h"
#include "FileIn.h"
#include "Rigid.h"
#include "LicenceSimple.h"

using Eigen::Vector3d;

int main(int argc, char *argv[]) {

	// タイトル表示．
	printf("        ---------------------------------------------------\n");
	printf("\n");
	printf("              DYNAMIC BEARING ANALYSIS IN NSK (D-BRAIN)\n");
	printf("\n");
	printf("                  D-B4P    : %s (Build Date)\n", __DATE__);
	printf("\n");
	printf("          Dynamic motion analysis of Ball 4 Point Bearing \n");
	printf("\n");
	printf("        ---------------------------------------------------\n");

	if (!LicenceSimple::hasTime(3000, 7, 30)) {
		std::cout << "エラー：ライセンスが切れています．" << std::endl;
		return -1;
	}

	char fpath_d4bin_csv[_MAX_STRING_];
	if (argc < 2) {
		std::cout << "インプットファイルの名前を入力してください" << std::endl;
		cin >> fpath_d4bin_csv;
	}
	else
		sprintf_s(fpath_d4bin_csv, _MAX_STRING_, argv[1]);

	B4P_FileIn FI;

	//CSVファイルから数値を読み込み，FIの各メンバ変数に書き込む
	bool has_error = FI.read_input_all(fpath_d4bin_csv);
	if(has_error) {
		std::cout <<  "inputファイル読み取り時にエラーが発生したため強制終了します．" << std::endl;
		return -1;
	}
	Rigid::l = FI.rigid.l;
	Rigid::t = FI.rigid.t;

	// このプログラムの指揮を司るオブジェクトの生成．
	B4P_Calculator calc;

	// 計算諸元を読み込んで初期化する．
	calc.b4p.init(FI, FI.DynSet.x, FI.DynSet.ax);

	// ファイル出力オブジェクトの作成．
	B4P_FileOut FO;
	FO.init(FI);

	// まずは静解析で力の釣り合いを解く．
	if (false) {

		// step 1. 各部材の位置の概算
		calc.Stt_init(FI.SttSet);
		double *x_stf = new double[calc.SttSet[0].n];
		calc.b4p.get_Xstf(x_stf);
		calc.Stt_solve_stf(x_stf);
		delete[] x_stf;
	}
	// 次に動解析で時間発展を計算する．
	calc.Dyn_init(FI.DynSet);
	double*y = new double[calc.dyn_n];
	double*dydt = new double[calc.dyn_n];
	double t[2];					// 動解析計算区間（時刻t[0]からt[1]まで計算）


	string Temp = "Temp.csv";
	// 続き計算モードではtempファイルから最新の時刻と座標を取得
	if (FI.DynSet.fromcontinuation) {
		double t_now;
		FO.read_Params(Temp, t_now, y);
		t[0] = t_now; t[1] = t_now + calc.t_step;
		calc.b4p.set_y(y);
	}
	// 初めから計算するモードではtempファイルの中身を全消去して新規作成
	else {
		FO.write_Header(Temp);
		t[0] = 0.; t[1] = calc.t_step;
		calc.b4p.get_y(y);
	}

	std::cout << "\n指定ループ回数を達成するか，Enterキーを押せば計算を終了します．" << std::endl << std::endl;

	std::chrono::system_clock::time_point start, end; // 型は auto で可
	start = std::chrono::system_clock::now(); // 計測開始時間

	if (FI.DynSet.stopcalculation) {
		// 1回目の計算で打ち切りしないようにトルクの初期値は十分大きい数値に設定
		double Ti0 = 1e10, Ti1 = 0, cnt = 0;
		while (t[0] < calc.calctime) {
			// odeintでyを時間発展させて解く．
			calc.Dyn_solve(y, t);

			// 次の時間ステップの設定．（t[0]はDyn_solve内で自動的にt[1]に変更される．）
			t[1] += calc.t_step;
			std::cout << "\r" << "Now " << 1e3 * t[0] << " [ms].\t"
				<< int(100.0 * t[0] / calc.calctime) << "% Complete." << std::string(20, ' ');

			// 1ステップ終了ごとにテキストに書き込む．
			calc.b4p.get_y(y);
			FO.write_Params(Temp, t[0], y);

			// enterキーを押したらループ中断
			if (_kbhit() != 0 && _getch() == '\r')
				break;

			calc.b4p.save(FO);
			Ti1 = FO.IR.T[0];
			double dTi = Unit::Nm2Nmm(Ti1 - Ti0);
			// トルクが非0かつ誤差閾値を下回ったらカウント
			if (abs(dTi) < FI.DynSet.dTierr && abs(Ti1) > 1e-9) {
				cnt++;
			}
			else {
				cnt = 0;
			}
			std::cout << "Ti = " << Unit::Nm2Nmm(Ti1) << "[Nmm]\tcount = " << cnt;
			if(cnt >= FI.DynSet.stp){
				std::cout << std::endl << std::endl << "打ち切り基準値を下回ったので計算を打ち切ります．" << std::endl;
				break;
			}
			Ti0 = Ti1;
		}
	}
	else {
		while (t[0] < calc.calctime) {
			// odeintでyを時間発展させて解く．
			calc.Dyn_solve(y, t);

			// 次の時間ステップの設定．（t[0]はDyn_solve内で自動的にt[1]に変更される．）
			t[1] += calc.t_step;
			std::cout << "\r" << "Now " << 1e3 * t[0] << " [ms].\t"
				<< int(100.0 * t[0] / calc.calctime) << "% Complete." << std::string(20, ' ');

			// 1ステップ終了ごとにテキストに書き込む．
			calc.b4p.get_y(y);
			FO.write_Params(Temp, t[0], y);

			// enterキーを押したらループ中断
			if (_kbhit() != 0 && _getch() == '\r')
				break;
		}
	}
	end = std::chrono::system_clock::now();  // 計測終了時間
	double elapsed = double(std::chrono::duration_cast<std::chrono::seconds>(end - start).count()); //処理に要した時間をミリ秒に変換

	std::cout << std::endl << std::endl << "計算に " << elapsed << " [s] かかりました．" << std::endl;
	std::cout << std::endl << "ポスト処理に移行します．" << std::endl;


	FO.write_AllHeader();

	// とりあえずベタ打ちでポスト処理．
	ifstream ifs(Temp, std::ios::in);
	string line;
	while (getline(ifs, line)) {
		vector<string> strvec = FileIn::split(line, ',');
		double*Params = new double[strvec.size()];
		for (int i = 0; i < int(strvec.size()) - 1; i++){
			Params[i] = stod(strvec.at(i + 1));
		}
		calc.b4p.set_y(Params);

		calc.b4p.get_dydt(y);
		calc.b4p.save(FO);
		FO.write_AllParams(stod(strvec.at(0)));
	}
	ifs.close();
	std::cout << std::endl << "全ての計算が完了しました．" << std::endl;

	// 動的確保した配列の解放．
	delete[] y, dydt;


	// 正常終了．
	return 0;

}


