#include "FileIn.h"


// csvファイルの各行をベクトルに読み込むためのメソッド．https://www.sejuku.net/blog/48660 からコピペ
vector<string> FileIn::split
(					// out:	入力文字列を区切り文字数(含む数:n)によってn+1個に分割した各文字列．
	string& input,	// in:	入力文字列．
	char delimiter	// in:	区切り文字．
) {
	istringstream stream(input);
	string field;
	vector<string> result;
	while (getline(stream, field, delimiter)) {
		result.push_back(field);
	}
	return result;
}

// csvファイルのうち$$で始まる行をstring型配列に格納するメソッド
vector<vector<string>>FileIn::input_to_array
(								// out:	入力ファイル全体の2次元配列．（大要素:各行．小要素:各列）
	const char filename[]		// in:	入力ファイル名．
) {
	//ファイルが開けなかったらエラー
	ifstream file(filename);
	if (!file) {
		cout << "Error! " << filename << " can not be opened" << endl;
		throw std::runtime_error("Inputfile Error");
	}

	// 各行に対して，(0)$$で始まれば2列目から読み取る(1)それ以外は読み飛ばす，のルールで(0)に該当したものだけ各配列に格納．
	string line;
	vector<vector<string>> inp_data;
	string dolmark = "$$";
	char delimiter = ',';
	while (getline(file, line)) {
		vector<string> strvec = FileIn::split(line, delimiter);
		size_t pos;

		//スペースを全部消去する(注意：以下のループを実施すると全角文字が文字化けを起こす)
		while ((pos = (strvec.at(0)).find_first_of("　 \t")) != string::npos) {
			strvec.at(0).erase(pos, 1);
		}

		//"$$"で始まる行であれば読み込む
		if (strvec.at(0).size() >= dolmark.size() && equal(begin(dolmark), end(dolmark), begin(strvec.at(0)))) {
			inp_data.push_back(strvec);
		}
	}
	return inp_data;
}


// 1列目に"param_name"が含まれている列を検索し，その行全体をstring配列に格納する．
vector<string> FileIn::pickup_data
(									// out:	検索で一番最初にヒットしたの行入力ファイル1次元配列．
	string param_name,				// in:	検索文字列．
	const vector<vector<string>>&inp_data	// in:	入力ファイル全体の2次元配列．（大要素:各行．小要素:各列）
) {
	for (unsigned int i = 0; i < inp_data.size(); ++i)
		if (inp_data[i][0] == param_name)
			return inp_data[i];

	//$$~~が見つからない時エラーを投げる
	throw std::runtime_error(param_name + " is Not Found");
}

// 1列目に"param_name"，2列目にelem_numが記載されている列を検索し，その行全体をstring配列に格納する．
vector<string> FileIn::pickup_multiple_data
(									// out:	検索で一番最初にヒットしたの行入力ファイル1次元配列．
	string param_name,				// in:	検索文字列．
	int elem_num,					// in:	要素番号．
	const vector<vector<string>>&inp_data	// in:	入力ファイル全体の2次元配列．（大要素:各行．小要素:各列）
) {
	for (unsigned int i = 0; i < inp_data.size(); ++i)
		if (inp_data[i][0]  == param_name &&
			stoi(inp_data[i][1]) == elem_num)
			return inp_data[i];

	//$$~~が見つからない時エラーを投げる
	throw std::runtime_error("NotFound");
}

// 1列目に"param_name"，2列目に"first_num"，3列目に"second_num"が記載されている列を検索し，その行全体をstring配列に格納する．
vector<string> FileIn::pickup_matrix_data
(									// out:	検索で一番最初にヒットしたの行入力ファイル1次元配列．
	string param_name,				// in:	検索文字列．
	int first_num,					// in:	要素番号．
	int second_num,					// in:	要素番号．
	const vector<vector<string>>&inp_data	// in:	入力ファイル全体の2次元配列．（大要素:各行．小要素:各列）
) {
	for (unsigned int i = 0; i < inp_data.size(); ++i)
		if (inp_data[i][0]  == param_name &&
			stoi(inp_data[i][1]) == first_num &&
			stoi(inp_data[i][2]) == second_num)
			return inp_data[i];

	//$$~~が見つからない時エラーを投げる
	throw std::runtime_error("NotFound");
}

// 入力ファイルの$$hoge列の行数に不足があった際に出すエラー．
bool FileIn::invalid_argument_error(const string & param_name) {
	cout <<  "エラー：：" << param_name << "の入力が誤っています．" << endl;
	cout <<  "　　　　　inputファイルを確認してください．" << endl;
	return true;
}

// 入力ファイルに$$hoge列が無かった際に出すエラー．
bool FileIn::cannot_find_error(const string & param_name) {
	cout <<  "エラー：：" << param_name << "の入力がありません．" << endl;
	cout <<  "　　　　　inputファイルを確認してください．" << endl;
	return true;
}



//while((pos = (strvec.at(0)).find_first_of(" ")) != string::npos){
			//if(strvec.at(0) == "$$"){
			//cout << strvec.at(0) << endl;
	//読み取りデータ確認用
	//for(unsigned int i = 0; i < inp_data.size(); i++){
	//	for(unsigned int j =0; j < inp_data.at(i).size();j++){
	//		cout << inp_data.at(i).at(j) << ", ";
	//	}
	//	cout << endl;
	//}
