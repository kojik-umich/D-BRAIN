/*******************************************************************************
!								"BS_FileIn.cpp"
!													2020/02/10	[Core-T]	楠崎
!
!*******************************************************************************/

#include "BS_FileIn.h"

void BS_FileIn::Nothing(void) {
}

//CSVファイルを読み込み，各メンバ変数に格納
void BS_FileIn::read_input_all(const char fpath_d4bin_csv[]) {

	//ファイルを読み込んでstring型の二次元ベクトルに書き込む
	vector<vector<string>> inp_data = FileIn::input_to_array(fpath_d4bin_csv);

	//ファイル読み取りエラーのフラッグ（true:エラーあり，false:エラー無し）
	bool has_error = false;

	// まず，動的確保に必要な情報を読み込み，メモリ確保まで行う．
	has_error |= this->read_allocate1(inp_data);
	has_error |= this->read_allocate2(inp_data);
	has_error |= this->read_allocate3(inp_data);

	has_error |= this->read_SttRotation(inp_data);

	// 現在の仕様ではとりあえず全ての玉径を一緒として扱う．
	has_error |= this->read_Ball(inp_data);

	// ナット読み込み．
	for (size_t i = 0; i < this->nut.size(); i++) {
		double Nut_PCD;
		has_error |= this->read_Nut(inp_data, i, Nut_PCD);
		for (size_t j = 0; j < this->nut[0].spiral.size(); j++)
			has_error |= this->read_NutSpiral(inp_data, i, j, Nut_PCD);
	}
	// 2つナットの場合は，与圧設定を読み込む．
	if (this->nut.size() == 2)
		has_error |= this->read_PreLoad(inp_data);
	this->read_BallNutPair(inp_data);

	// シャフト読み込み．
	double Shaft_PCD;
	has_error |= this->read_Shaft(inp_data, Shaft_PCD);
	for (size_t i = 0; i < this->shaft[0].spiral.size(); i++)
		has_error |= this->read_ShaftSpiral(inp_data, i, Shaft_PCD);
	this->read_BallShaftPair(inp_data);
	
	// 位相角の計算（$$Circuitの読み込みの時点では計算できないため，ここで計算）
	for (int i = 0; i < this->circuit.size(); i++) 
		this->set_Circuitphase(i);
	this->read_ShaftMassSet(inp_data);
	// 荷重読み込み．
	for (size_t i = 0; i < this->load.size(); i++)
		has_error |= this->read_SttLoad(inp_data, i);

	// 以下，各パラメタに値を代入していく．
	has_error |= this->read_RollingResistance(inp_data);
	has_error |= this->read_Coulomb(inp_data);
	has_error |= this->read_Ellipse(inp_data);
	has_error |= this->read_FilmThickness(inp_data);
	has_error |= this->read_Dimension(inp_data);
	has_error |= this->read_Oil(inp_data);
	has_error |= this->read_Gravity(inp_data);
	has_error |= this->read_Bound(inp_data);
	has_error |= this->read_PositionSet(inp_data);
	if (this->initial.preset == Initial::ReadPos)
		has_error |= this->read_Position(inp_data);

	has_error |= this->read_SttMode(inp_data);
	for (int i = 0; i < 3; i++)
		has_error |= this->read_SttSet(inp_data, i);

	for (int i = 0; i < 2; i++)
		has_error |= this->read_DynSet(inp_data, i);
	has_error |= this->read_Output(inp_data);


	has_error |= this->read_DynRotationStep(inp_data);
	for (size_t i = 0; i < this->dyn.wxt_n; i++)
		has_error |= this->read_DynRotation(inp_data, i);

	has_error |= this->read_DynVelocityStep(inp_data);
	for (size_t i = 0; i < this->dyn.vxd_n; i++)
		has_error |= this->read_DynVelocity(inp_data, i);

	has_error |= this->read_DynLoadStep(inp_data);
	for (size_t i = 0; i < this->dyn.ft_n; i++)
		has_error |= this->read_DynLoad(inp_data, i);

	for (int i = 0; i < 2; i++)
		has_error |= this->read_StopDynCalc(inp_data, i);

	has_error |= this->read_DeleteLastOutput(inp_data);

	
	// エラーを検知したら，例外をスローする．
	if (has_error)
		throw std::runtime_error("Inputfile Error");

	return;
}

// 動的確保（new）する部分だけ読み込み・実行するサブルーチン．その１．
bool BS_FileIn::read_allocate1(const vector<vector<string>>&inp_data) {

	bool has_error = false;

	this->shaft.resize(1);

	int spiralnum;
	has_error |= this->read_SpiralNum(inp_data, spiralnum);
	this->shaft[0].spiral.resize(spiralnum);

	int nutnum;
	has_error |= this->read_NutNum(inp_data, nutnum);
	this->nut.resize(nutnum);
	for (int i = 0; i < nutnum; i++)
		this->nut[i].spiral.resize(spiralnum);

	int circuitnum;
	has_error |= this->read_CircuitNum(inp_data, circuitnum);
	this->circuit.resize(circuitnum);

	return has_error;
}

// 動的確保（new）する部分だけ読み込み・実行するサブルーチン．その２．
bool BS_FileIn::read_allocate2(const vector<vector<string>>&inp_data) {

	bool has_error = false;

	this->ballnum = 0;
	for (size_t i = 0; i < this->circuit.size(); i++) {
		int ballnum;
		has_error |= this->read_Circuit(inp_data, i, ballnum);
		this->circuit[i].ball.resize(ballnum);
		this->ballnum += ballnum;
	}
	this->BallNutPair.resize(this->ballnum);
	this->BallShaftPair.resize(this->ballnum);

	return has_error;
}

// 動的確保（new）する部分だけ読み込み・実行するサブルーチン．その３．
bool BS_FileIn::read_allocate3(const vector<vector<string>>&inp_data) {

	bool has_error = false;
	int sttloadnum;
	has_error |= this->read_SttLoadNum(inp_data, sttloadnum);
	this->load.resize(sttloadnum);

	return has_error;
}

// ナット数の読み込み．
bool BS_FileIn::read_NutNum(const vector<vector<string>>&inp_data, int&nutnum) {
	string param_name = "$$NutNum";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		nutnum = stoi(param[1]);
		if (nutnum < 1 || nutnum > 2) {
			cout << "ナット数の指定が" << nutnum << "になっています．" << endl;
			cout << "ナット数は1か2を入力してください．" << endl;
			return true;
		}
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// ねじ条数の読み込み．
bool BS_FileIn::read_SpiralNum(const vector<vector<string>>&inp_data, int&spiralnum) {
	string param_name = "$$SpiralNum";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		spiralnum = stoi(param[1]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 回路数の読み込み．
bool BS_FileIn::read_CircuitNum(const vector<vector<string>>&inp_data, int&circuitnum) {
	string param_name = "$$CircuitNum";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		circuitnum = stoi(param[1]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 回路情報の読み込み．
bool BS_FileIn::read_Circuit(const vector<vector<string>>&inp_data, int i, int&ballnum) {
	string param_name = "$$Circuit";
	try {
		vector<string> param = FileIn::pickup_multiple_data(param_name, i, inp_data);
		this->circuit[i].inut = stoi(param[2]);
		this->circuit[i].is = stoi(param[3]);
		ballnum = stoi(param[4]);
		// 負荷開始点・終了点はあとで位相角に変換
		this->circuit[i].x0 =  Unit::mm2m(stod(param[5]));
		this->circuit[i].x1 =  Unit::mm2m(stod(param[6]));
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 回路の位相角を計算し，構造体に代入
void BS_FileIn::set_Circuitphase(int i) {

	int inut = this->circuit[i].inut;
	int is  = this->circuit[i].is;
	this->circuit[i].th0 = calc_phase(this->nut[inut].spiral[is].l, this->circuit[i].x0);
	this->circuit[i].th1 = calc_phase(this->nut[inut].spiral[is].l, this->circuit[i].x1);
	return;

}

// x座標を位相角に変換
double BS_FileIn::calc_phase(
	double l,	// リード長[m]
	double x	// x座標[m]
) {
	double phase = x / l * 2 * Numeric::pi;
	return phase;
}

// 与圧方式の読み込み．
bool BS_FileIn::read_PreLoad(const vector<vector<string>>&inp_data) {
	string param_name = "$$PreLoad";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->preload.mode = static_cast<Preload::Mode>(stoi(param[1]));
		switch (this->preload.mode) {
		case Preload::Mode::distance:
			this->nut[0].x0 = Unit::mm2m(stod(param[2]));
			this->nut[1].x0 = Unit::mm2m(stod(param[3]));
			break;
		case Preload::Mode::load:
			cout << "まだ定圧与圧モードは実装していません．" << endl;
			return true;
		default:
			cout << "与圧方式は0か1で入力してください．" << endl;
			return true;
		}
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 転がり摩擦設定条件データ取得
bool BS_FileIn::read_RollingResistance(const vector<vector<string>>&inp_data) {
	string param_name = "$$RollingResistance";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->tribology.rollingresistance = static_cast<BS_In::Tribology::RollingResistance>(stoi(param[1]));
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 金属接触滑り摩擦設定条件データ取得
bool BS_FileIn::read_Coulomb(const vector<vector<string>>&inp_data) {
	string param_name = "$$Coulomb";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->tribology.coulomb = static_cast<BS_In::Tribology::Coulomb>(stoi(param[1]));
		this->tribology.coulomb_slope = stod(param[2]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 金属接触滑り摩擦設定条件データ取得
bool BS_FileIn::read_Ellipse(const vector<vector<string>>&inp_data) {
	string param_name = "$$Ellipse";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->tribology.ellipse= static_cast<BS_In::Tribology::Ellipse>(stoi(param[1]));
		this->tribology.ellipse_mesh[0] = stoi(param[2]);
		this->tribology.ellipse_mesh[1] = stoi(param[3]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 油膜厚さ計算式データ取得
bool BS_FileIn::read_FilmThickness(const vector<vector<string>>&inp_data) {
	string param_name = "$$FilmThickness";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->tribology.filmThickness = static_cast<BS_In::Tribology::FilmThickness>(stoi(param[1]));
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// ヒステリシス損失データ取得
bool BS_FileIn::read_Hysteresis(const vector<vector<string>>&inp_data) {
	string param_name = "$$Hysteresis";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->tribology.hysteresis = static_cast<Tribology::Hysteresis>(stoi(param[1]));
		this->tribology.hysteresis_factor = stod(param[2]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

bool BS_FileIn::read_Dimension(const vector<vector<string>>&inp_data) {
	string param_name = "$$Dimension";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->rigid.l = stod(param[1]);
		this->rigid.t = stod(param[2]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// ボール入力．
bool BS_FileIn::read_Ball(const vector<vector<string>>&inp_data) {
	string param_name = "$$Ball";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		for (size_t i = 0; i < this->circuit.size(); i++) {
			for (size_t j = 0; j < this->circuit[i].ball.size(); j++) {
				this->circuit[i].ball[j].density = stod(param[1]);
				this->circuit[i].ball[j].young = Unit::GPa2Pa(stod(param[2]));
				this->circuit[i].ball[j].poisson = stod(param[3]);
				this->circuit[i].ball[j].rms = Unit::um2m(stod(param[4]));
				this->circuit[i].ball[j].r = Unit::mm2m(stod(param[5])) / 2;
			}
		}
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// シャフト入力．
bool BS_FileIn::read_Shaft(const vector<vector<string>>&inp_data, double&Shaft_PCD) {
	string param_name = "$$Shaft";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->shaft[0].density = stod(param[1]);
		this->shaft[0].young = Unit::GPa2Pa(stod(param[2]));
		this->shaft[0].poisson = stod(param[3]);
		Shaft_PCD = Unit::mm2m(stod(param[4]));
		this->shaft[0].ri = Unit::mm2m(stod(param[5])) / 2;
		this->shaft[0].ro = Unit::mm2m(stod(param[6])) / 2;
		this->shaft[0].l = Unit::mm2m(stod(param[7]));

		if (param[8] == "auto")
			this->shaft[0].m = BS_FileIn::calc_Cylinder_m(this->shaft[0].ri, this->shaft[0].ro, this->shaft[0].l, this->shaft[0].density);
		else
			this->shaft[0].m = stod(param[8]);
		if (param[9] == "auto")
			this->shaft[0].Ix = BS_FileIn::calc_Cylinder_Ix(this->shaft[0].ri, this->shaft[0].ro, this->shaft[0].l, this->shaft[0].density);
		else
			this->shaft[0].Ix = stod(param[9]);
		if (param[10] == "auto")
			this->shaft[0].Iyz = BS_FileIn::calc_Cylinder_Iyz(this->shaft[0].ri, this->shaft[0].ro, this->shaft[0].l, this->shaft[0].density);
		else
			this->shaft[0].Iyz = stod(param[10]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// シャフト螺旋入力．
bool BS_FileIn::read_ShaftSpiral(const vector<vector<string>>&inp_data, int i, double PCD) {
	string param_name = "$$ShaftSpiral";
	try {
		vector<string> param = FileIn::pickup_multiple_data(param_name, i, inp_data);
		this->shaft[0].spiral[i].alp = Unit::deg2rad(stod(param[2]));
		this->shaft[0].spiral[i].l = Unit::mm2m(stod(param[3]));
		this->shaft[0].spiral[i].r = PCD / 2;

		for (int j = 0; j < 2; j++) {
			int k = 4 * j + 4;
			this->shaft[0].spiral[i].groove[j].sigma = Unit::um2m(stod(param[k + 0]));
			this->shaft[0].spiral[i].groove[j].r = Unit::mm2m(stod(param[k + 1]));
			this->shaft[0].spiral[i].groove[j].eta[0] = Unit::mm2m(stod(param[k + 2]));
			this->shaft[0].spiral[i].groove[j].eta[1] = Unit::mm2m(stod(param[k + 3]));

		}
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

bool BS_FileIn::read_BallShaftPair(const vector<vector<string>>&inp_data) {

	string param_name = "$$BallShaftPair";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);

		for (size_t i = 0; i < this->BallShaftPair.size(); i++)
			for (size_t j = 0; j < 2; j++) {
				this->BallShaftPair[i].groove[j].mu = stod(param[2 * j + 1]);
				this->BallShaftPair[i].groove[j].zeta = stod(param[2 * j + 2]);
			}
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// ナット入力．
bool BS_FileIn::read_Nut(const vector<vector<string>>&inp_data, int i, double&Nut_PCD) {
	string param_name = "$$Nut";
	try {
		vector<string> param = FileIn::pickup_multiple_data(param_name, i, inp_data);
		this->nut[i].density = stod(param[2]);
		this->nut[i].young = Unit::GPa2Pa(stod(param[3]));
		this->nut[i].poisson = stod(param[4]);
		Nut_PCD = Unit::mm2m(stod(param[5]));
		this->nut[i].ri = Unit::mm2m(stod(param[6])) / 2;
		this->nut[i].ro = Unit::mm2m(stod(param[7])) / 2;
		this->nut[i].l = Unit::mm2m(stod(param[8]));
		this->nut[i].xg = Unit::mm2m(stod(param[9]));

		if (param[10] == "auto")
			this->nut[i].m = BS_FileIn::calc_Cylinder_m(this->nut[i].ri, this->nut[i].ro, this->nut[i].l, this->nut[i].density);
		else
			this->nut[i].m = stod(param[10]);

		if (param[11] == "auto")
			this->nut[i].Ix = BS_FileIn::calc_Cylinder_Ix(this->nut[i].ri, this->nut[i].ro, this->nut[i].l, this->nut[i].density);
		else
			this->nut[i].Ix = stod(param[11]);

		if (param[12] == "auto")
			this->nut[i].Iyz = BS_FileIn::calc_Cylinder_Iyz(this->nut[i].ri, this->nut[i].ro, this->nut[i].l, this->nut[i].density);
		else
			this->nut[i].Iyz = stod(param[12]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

double BS_FileIn::calc_Cylinder_m(double ri, double ro, double l, double rho) {
	double Si = Numeric::Square(ri) * Numeric::pi;
	double So = Numeric::Square(ro) * Numeric::pi;
	double m = (So - Si) * l * rho;
	return m;
}

double BS_FileIn::calc_Cylinder_Ix(double ri, double ro, double l, double rho) {
	double m = BS_FileIn::calc_Cylinder_m(ri, ro, l, rho);
	double Di2 = Numeric::Square(ri * 2);
	double Do2 = Numeric::Square(ro * 2);
	double Ix = 1.0 / 8 * m * (Di2 + Do2);
	return Ix;
}

double BS_FileIn::calc_Cylinder_Iyz(double ri, double ro, double l, double rho) {
	double m = BS_FileIn::calc_Cylinder_m(ri, ro, l, rho);
	double Di2 = Numeric::Square(ri * 2);
	double Do2 = Numeric::Square(ro * 2);
	double Iyz = 1.0 / 4 * m * ((Di2 + Do2) / 4 + Numeric::Square(l) / 3);
	return Iyz;
}

double BS_FileIn::calc_Cylinder_Iyz_l(double m, double x, double l) {
	double Iyz = m * (3 * x * x + l * l / 4);
	return Iyz;
}

// ナット螺旋入力．
bool BS_FileIn::read_NutSpiral(const vector<vector<string>>&inp_data, int i, int j, double PCD) {

	string param_name = "$$NutSpiral";
	try {
		vector<string> param = FileIn::pickup_matrix_data(param_name, i, j, inp_data);
		this->nut[i].spiral[j].alp = Unit::deg2rad(stod(param[3]));
		this->nut[i].spiral[j].l = Unit::mm2m(stod(param[4]));
		this->nut[i].spiral[j].r = PCD / 2;
		for (int k = 0; k < 2; k++) {
			int l = 4 * k + 5;
			this->nut[i].spiral[j].groove[k].sigma = Unit::um2m(stod(param[l + 0]));
			this->nut[i].spiral[j].groove[k].r = Unit::mm2m(stod(param[l + 1]));
			this->nut[i].spiral[j].groove[k].eta[0] = Unit::mm2m(stod(param[l + 2]));
			this->nut[i].spiral[j].groove[k].eta[1] = Unit::mm2m(stod(param[l + 3]));
		}
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

bool BS_FileIn::read_BallNutPair(const vector<vector<string>>&inp_data) {

	string param_name = "$$BallNutPair";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);

		for (size_t i = 0; i < this->BallNutPair.size(); i++)
			for (size_t j = 0; j < 2; j++) {
				this->BallNutPair[i].groove[j].mu = stod(param[2 * j + 1]);
				this->BallNutPair[i].groove[j].zeta = stod(param[2 * j + 2]);
			}
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}


// 油物性入力．
bool BS_FileIn::read_Oil(const vector<vector<string>>&inp_data) {
	string param_name = "$$Oil";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->oil.alpha = Unit::mm2kgf2Painv(stod(param[1]));
		this->oil.beta = stod(param[2]);
		this->oil.k = stod(param[3]);
		this->oil.eta = stod(param[4]);
		this->oil.lm = stod(param[5]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 重力入力．
bool BS_FileIn::read_Gravity(const vector<vector<string>>&inp_data) {
	string param_name = "$$Gravity";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->rigid.g[0] = stod(param[1]);
		this->rigid.g[1] = stod(param[2]);
		this->rigid.g[2] = stod(param[3]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// シャフト位置入力．
bool BS_FileIn::read_PositionSet(const vector<vector<string>>&inp_data) {
	string param_name = "$$PositionSet";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->initial.preset = static_cast<BS_FileIn::Initial::Preset>(stoi(param[1]));
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// シャフト初期値or拘束値入力．
bool BS_FileIn::read_Position(const vector<vector<string>>&inp_data) {
	string param_name = "$$Position";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->bound.x0[0] = Unit::mm2m(stod(param[1]));
		this->bound.x0[1] = Unit::mm2m(stod(param[2]));
		this->bound.x0[2] = Unit::mm2m(stod(param[3]));
		double thy = Unit::deg2rad(stod(param[4]));
		double thz = Unit::deg2rad(stod(param[5]));

		this->bound.ax0[0] = 1.0;
		this->bound.ax0[1] = tan(thz);
		this->bound.ax0[2] =-tan(thy);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 荷重数の読み込み．
bool BS_FileIn::read_SttLoadNum(const vector<vector<string>>&inp_data, int&sttloadnum) {
	string param_name = "$$SttLoadNum";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		sttloadnum = stoi(param[1]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 荷重入力．
bool BS_FileIn::read_SttLoad(const vector<vector<string>>&inp_data, int i) {
	string param_name = "$$SttLoad";
	try {
		vector<string> param = FileIn::pickup_multiple_data(param_name, i, inp_data);
		this->load[i].x[0] = Unit::mm2m(stod(param[2]));
		this->load[i].x[1] = Unit::mm2m(stod(param[3]));
		this->load[i].x[2] = Unit::mm2m(stod(param[4]));
		this->load[i].F[0] = stod(param[5]);
		this->load[i].F[1] = stod(param[6]);
		this->load[i].F[2] = stod(param[7]);
		this->load[i].T[0] = Unit::Nmm2Nm(stod(param[8]));
		this->load[i].T[1] = Unit::Nmm2Nm(stod(param[9]));
		this->load[i].T[2] = Unit::Nmm2Nm(stod(param[10]));
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 回転数入力．
bool BS_FileIn::read_SttRotation(const vector<vector<string>>&inp_data) {
	string param_name = "$$SttRotation";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->stt.w0 = Unit::rpm2rads(stod(param[1]));
		this->stt.wn = Unit::rpm2rads(stod(param[2]));
		double l0 = Unit::mm2m(stod(param[3]));
		this->stt.v0 = (this->stt.w0 - this->stt.wn) * l0 / 2 / Numeric::pi;
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 静解析の釣り合い条件．
bool BS_FileIn::read_SttMode(const vector<vector<string>>&inp_data) {
	string param_name = "$$SttMode";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->runStatic  = stoi(param[1]) != 0;
		this->stt.run[0] = stoi(param[2]) != 0;
		this->stt.run[1] = stoi(param[3]) != 0;
		this->stt.run[2] = stoi(param[4]) != 0;
		this->stt.run[3] = stoi(param[5]) != 0;
		this->stt.run[4] = stoi(param[6]) != 0;
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 静解析条件入力．
bool BS_FileIn::read_SttSet(const vector<vector<string>>&inp_data, int i) {
	string param_name = "$$SttSet" + to_string(i);
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->stt.set[i].iter1 = stoi(param[1]);
		this->stt.set[i].iter2 = stoi(param[2]);
		this->stt.set[i].rs = stod(param[3]);
		this->stt.set[i].jac_eps = stod(param[4]);
		for (int j = 0; j < 6; j++)
			this->stt.set[i].eps[j] = stod(param[j + 5]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 動解析条件入力．
bool BS_FileIn::read_DynSet(const vector<vector<string>>&inp_data, int i) {
	string param_name = "$$DynSet" + to_string(i);
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->dyn.set[i].calctime = Unit::ms2s(stod(param[1]));
		this->dyn.set[i].sampling = Numeric::Roundoff(stod(param[2]));
		this->dyn.set[i].h = stod(param[3]);
		this->dyn.set[i].hmin = stod(param[4]);
		this->dyn.set[i].ep = stod(param[5]);
		this->dyn.set[i].tr = stod(param[6]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

bool BS_FileIn::read_Bound(const vector<vector<string>>& inp_data) {
	string param_name = "$$Bound";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->bound.v_const[0] = stoi(param[1]) != 0;
		this->bound.v_const[1] = stoi(param[2]) != 0;
		this->bound.v_const[2] = stoi(param[3]) != 0;
		this->bound.w_const[0] = stoi(param[4]) != 0;
		this->bound.w_const[1] = stoi(param[5]) != 0;
		this->bound.w_const[2] = stoi(param[6]) != 0;
	}					   
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// $$ShaftMassSet 読み込み
bool BS_FileIn::read_ShaftMassSet(const vector<vector<string>>& inp_data) {
	string param_name = "$$ShaftMassSet";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->mass = static_cast<BS_FileIn::Mass>(stoi(param[1]));
		this->calc_Mass();
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}
void BS_FileIn::calc_Mass(void) {
	using namespace Tribology;
	double nut_m = 0.0;
	double nut_Ix = 0.0;
	double nut_Iyz = 0.0;
	for (size_t i = 0; i < this->nut.size(); i++) {
		nut_m += this->nut[i].m;
		nut_Ix += this->nut[i].Ix;
		nut_Iyz += BS_FileIn::calc_Cylinder_Iyz_l(this->nut[i].m, this->nut[i].xg, this->nut[i].l);
	}
	switch (this->mass) {
	case BS_FileIn::Mass::Reduced:
		this->shaft[0].m = ReducedMass(nut_m, this->shaft[0].m);
		this->shaft[0].Ix = ReducedMass(nut_Ix, this->shaft[0].Ix);
		this->shaft[0].Iyz = ReducedMass(nut_Iyz, this->shaft[0].Iyz);
		break;
	case BS_FileIn::Mass::Nut:
		this->shaft[0].m = nut_m;
		this->shaft[0].Ix = nut_Ix;
		this->shaft[0].Iyz = nut_Iyz;
		break;
	case BS_FileIn::Mass::Shaft:
		this->shaft[0].m = this->shaft[0].m;
		this->shaft[0].Ix = this->shaft[0].Ix;
		this->shaft[0].Iyz = this->shaft[0].Iyz;
		break;
	}
	return;
}

// 出力ファイル名入力．
bool BS_FileIn::read_Output(const vector<vector<string>>&inp_data) {
	string param_name = "$$Output";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->output.temp = param[1];
		this->output._01 = param[2];
		this->output._02 = param[3];
		this->output._03 = param[4];
		this->output._04 = param[5];
		this->output._05 = param[6];
		this->output._06 = param[7];
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 動解析詳細入力条件読み取り．
bool BS_FileIn::read_DynRotationStep(const vector<vector<string>>&inp_data) {
	string param_name = "$$DynRotationStep";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->dyn.wxt_n = stoi(param[1]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 動解析詳細入力条件読み取り．
bool BS_FileIn::read_DynRotation(const vector<vector<string>>&inp_data, int i) {
	string param_name = "$$DynRotation";
	try {
		vector<string> param = FileIn::pickup_multiple_data(param_name, i, inp_data);
		double t = Unit::ms2s(stod(param[2]));
		this->dyn.wxt[t] = Unit::rpm2rads(stod(param[3]));
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}


// 動解析詳細入力条件読み取り．
bool BS_FileIn::read_DynVelocityStep(const vector<vector<string>>&inp_data) {
	string param_name = "$$DynVelocityStep";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->dyn.vxd_n = stoi(param[1]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 動解析詳細入力条件読み取り．
bool BS_FileIn::read_DynVelocity(const vector<vector<string>>&inp_data, int i) {
	string param_name = "$$DynVelocity";
	try {
		vector<string> param = FileIn::pickup_multiple_data(param_name, i, inp_data);
		double t = Unit::ms2s(stod(param[2]));
		this->dyn.vxt[t] = Unit::rpm2rads(stod(param[3]));
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}


// 動解析詳細入力条件読み取り．
bool BS_FileIn::read_DynLoadStep(const vector<vector<string>>&inp_data) {
	string param_name = "$$DynLoadStep";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->dyn.ft_n = stoi(param[1]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 動解析詳細入力条件読み取り．
bool BS_FileIn::read_DynLoad(const vector<vector<string>>&inp_data, int i) {
	string param_name = "$$DynLoad";
	try {
		vector<string> param = FileIn::pickup_multiple_data(param_name, i, inp_data);
		double t = Unit::ms2s(stod(param[2]));
		vector<double> Ft(6);
		Ft[0] = stod(param[3]);
		Ft[1] = stod(param[4]);
		Ft[2] = stod(param[5]);
		Ft[3] = stod(param[6]);
		Ft[4] = stod(param[7]);
		Ft[5] = stod(param[8]);
		this->dyn.ft[t] = Ft;
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 動解析計算打ち切り設定
bool BS_FileIn::read_StopDynCalc(const vector<vector<string>>&inp_data, int i) {
	string param_name = "$$StopDynCalc" + to_string(i);
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->dyn.set[i].stopcalc = (stoi(param[1]) != 0);
		this->dyn.set[i].dTerr = Unit::Nmm2Nm(stod(param[2]));
		this->dyn.set[i].stp = stoi(param[3]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// 前回計算結果消去
bool BS_FileIn::read_DeleteLastOutput(const vector<vector<string>>&inp_data) {
	string param_name = "$$DeleteLastOutput";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->output.deletelastout = (stoi(param[1]) != 0);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}



