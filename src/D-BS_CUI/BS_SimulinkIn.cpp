/*******************************************************************************
!								"BS_SimulinkIn.cpp"
!													2020/02/10	[Core-T]	楠崎
!
!*******************************************************************************/

#include "BS_SimulinkIn.h"

void BS_SimulinkIn::Nothing(void) {
}

//CSVファイルを読み込み，各メンバ変数に格納
void BS_SimulinkIn::read_csv(const char fpath_d4bin_csv[]) {

	//ファイルを読み込んでstring型の二次元ベクトルに書き込む
	vector<vector<string>> inp_data = FileIn::input_to_array(fpath_d4bin_csv);

	//ファイル読み取りエラーのフラッグ（true:エラーあり，false:エラー無し）
	bool has_error = false;

	// まず，動的確保に必要な情報を読み込み，メモリ確保まで行う．
	has_error |= this->read_allocate1(inp_data);
	has_error |= this->read_allocate2(inp_data);

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

	// 以下，各パラメタに値を代入していく．
	has_error |= this->read_RollingResistance(inp_data);
	has_error |= this->read_Coulomb(inp_data);
	has_error |= this->read_FilmThickness(inp_data);
	has_error |= this->read_Dimension(inp_data);
	has_error |= this->read_Oil(inp_data);
	has_error |= this->read_Gravity(inp_data);

	// エラーを検知したら，例外をスローする．
	if (has_error)
		throw std::runtime_error("Inputfile Error");

	return;
}

// 動的確保（new）する部分だけ読み込み・実行するサブルーチン．その１．
bool BS_SimulinkIn::read_allocate1(const vector<vector<string>>&inp_data) {

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
bool BS_SimulinkIn::read_allocate2(const vector<vector<string>>&inp_data) {

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

// ナット数の読み込み．
bool BS_SimulinkIn::read_NutNum(const vector<vector<string>>&inp_data, int&nutnum) {
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
bool BS_SimulinkIn::read_SpiralNum(const vector<vector<string>>&inp_data, int&spiralnum) {
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
bool BS_SimulinkIn::read_CircuitNum(const vector<vector<string>>&inp_data, int&circuitnum) {
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
bool BS_SimulinkIn::read_Circuit(const vector<vector<string>>&inp_data, int i, int&ballnum) {
	string param_name = "$$Circuit";
	try {
		vector<string> param = FileIn::pickup_multiple_data(param_name, i, inp_data);
		this->circuit[i].inut = stoi(param[2]);
		this->circuit[i].is = stoi(param[3]);
		ballnum = stoi(param[4]);
		// 負荷開始点・終了点はあとで位相角に変換
		this->circuit[i].x0 = Unit::mm2m(stod(param[5]));
		this->circuit[i].x1 = Unit::mm2m(stod(param[6]));
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
void BS_SimulinkIn::set_Circuitphase(int i) {

	int inut = this->circuit[i].inut;
	int is = this->circuit[i].is;
	this->circuit[i].th0 = calc_phase(this->nut[inut].spiral[is].l, this->circuit[i].x0);
	this->circuit[i].th1 = calc_phase(this->nut[inut].spiral[is].l, this->circuit[i].x1);
	return;

}

// x座標を位相角に変換
double BS_SimulinkIn::calc_phase(
	double l,	// リード長[m]
	double x	// x座標[m]
) {
	double phase = x / l * 2 * Numeric::pi;
	return phase;
}

// 与圧方式の読み込み．
bool BS_SimulinkIn::read_PreLoad(const vector<vector<string>>&inp_data) {
	string param_name = "$$PreLoad";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->nut[0].x0 = Unit::mm2m(stod(param[2]));
		this->nut[1].x0 = Unit::mm2m(stod(param[3]));
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
bool BS_SimulinkIn::read_RollingResistance(const vector<vector<string>>&inp_data) {
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
bool BS_SimulinkIn::read_Coulomb(const vector<vector<string>>&inp_data) {
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

// 油膜厚さ計算式データ取得
bool BS_SimulinkIn::read_FilmThickness(const vector<vector<string>>&inp_data) {
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
bool BS_SimulinkIn::read_Hysteresis(const vector<vector<string>>&inp_data) {
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

bool BS_SimulinkIn::read_Dimension(const vector<vector<string>>&inp_data) {
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
bool BS_SimulinkIn::read_Ball(const vector<vector<string>>&inp_data) {
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
bool BS_SimulinkIn::read_Shaft(const vector<vector<string>>&inp_data, double&Shaft_PCD) {
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
			this->shaft[0].m = BS_SimulinkIn::calc_Cylinder_m(this->shaft[0].ri, this->shaft[0].ro, this->shaft[0].l, this->shaft[0].density);
		else
			this->shaft[0].m = stod(param[8]);
		if (param[9] == "auto")
			this->shaft[0].Ix = BS_SimulinkIn::calc_Cylinder_Ix(this->shaft[0].ri, this->shaft[0].ro, this->shaft[0].l, this->shaft[0].density);
		else
			this->shaft[0].Ix = stod(param[9]);
		if (param[10] == "auto")
			this->shaft[0].Iyz = BS_SimulinkIn::calc_Cylinder_Iyz(this->shaft[0].ri, this->shaft[0].ro, this->shaft[0].l, this->shaft[0].density);
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
bool BS_SimulinkIn::read_ShaftSpiral(const vector<vector<string>>&inp_data, int i, double PCD) {
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

bool BS_SimulinkIn::read_BallShaftPair(const vector<vector<string>>&inp_data) {

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
bool BS_SimulinkIn::read_Nut(const vector<vector<string>>&inp_data, int i, double&Nut_PCD) {
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
			this->nut[i].m = BS_SimulinkIn::calc_Cylinder_m(this->nut[i].ri, this->nut[i].ro, this->nut[i].l, this->nut[i].density);
		else
			this->nut[i].m = stod(param[10]);

		if (param[11] == "auto")
			this->nut[i].Ix = BS_SimulinkIn::calc_Cylinder_Ix(this->nut[i].ri, this->nut[i].ro, this->nut[i].l, this->nut[i].density);
		else
			this->nut[i].Ix = stod(param[11]);

		if (param[12] == "auto")
			this->nut[i].Iyz = BS_SimulinkIn::calc_Cylinder_Iyz(this->nut[i].ri, this->nut[i].ro, this->nut[i].l, this->nut[i].density);
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

double BS_SimulinkIn::calc_Cylinder_m(double ri, double ro, double l, double rho) {
	double Si = Numeric::Square(ri) * Numeric::pi;
	double So = Numeric::Square(ro) * Numeric::pi;
	double m = (So - Si) * l * rho;
	return m;
}

double BS_SimulinkIn::calc_Cylinder_Ix(double ri, double ro, double l, double rho) {
	double m = BS_SimulinkIn::calc_Cylinder_m(ri, ro, l, rho);
	double Di2 = Numeric::Square(ri * 2);
	double Do2 = Numeric::Square(ro * 2);
	double Ix = 1.0 / 8 * m * (Di2 + Do2);
	return Ix;
}

double BS_SimulinkIn::calc_Cylinder_Iyz(double ri, double ro, double l, double rho) {
	double m = BS_SimulinkIn::calc_Cylinder_m(ri, ro, l, rho);
	double Di2 = Numeric::Square(ri * 2);
	double Do2 = Numeric::Square(ro * 2);
	double Iyz = 1.0 / 4 * m * ((Di2 + Do2) / 4 + Numeric::Square(l) / 3);
	return Iyz;
}

double BS_SimulinkIn::calc_Cylinder_Iyz_l(double m, double x, double l) {
	double Iyz = m * (3 * x * x + l * l / 4);
	return Iyz;
}

// ナット螺旋入力．
bool BS_SimulinkIn::read_NutSpiral(const vector<vector<string>>&inp_data, int i, int j, double PCD) {

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

bool BS_SimulinkIn::read_BallNutPair(const vector<vector<string>>&inp_data) {

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
bool BS_SimulinkIn::read_Oil(const vector<vector<string>>&inp_data) {
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
bool BS_SimulinkIn::read_Gravity(const vector<vector<string>>&inp_data) {
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

