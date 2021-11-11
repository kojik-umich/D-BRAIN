#include "FileOut.h"


// ファイルを新規作成し，ヘッダーを書き込むメソッド．
void FileOut::write_header(
	string name,	// [in]	ファイル名
	string head		// [in]	書き込むヘッダー文字列
) {
	ofstream file;
	file.open(name.c_str(), ios_base::trunc);
	file << head;

	return;
}

// 作成されているファイルに追加で値を1行に書き込むメソッド．
void FileOut::write_vector(
	string name,			// [in]	ファイル名
	const VectorXd&vec,		// [in]	書き込む数値配列
	int prec				// [in]	数値の有効数字
) {
	IOFormat CSVFormat(prec, 0, ", ", ", ", "", "", "", "\n");
	ofstream file;
	file.open(name.c_str(), ios_base::app);
	file << vec.format(CSVFormat);

	return;
}
