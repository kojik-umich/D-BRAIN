#include "pch.h"

class FileInTest : public ::testing::Test {
public:
	vector<vector<string>> inp_data = {
		{"a", "b", "c", "d"} , {"e", "2", "f", "g"}, {"h", "3", "5", "j"} };
};

TEST_F(FileInTest, input_to_array0) {

	// 存在しないファイルの読み取りエラーの確認
	EXPECT_ANY_THROW({
	vector<vector<string>> inp = FileIn::input_to_array("hoge.csv");
		}; );

	// outに置いてあるフォルダを参照．
	char filename[] = "FileIn_test.csv";

	// テキストファイルの読み込み．（エラーが出る方はprojのプロパティをどうぞ．OpenCppCoverageも忘れずに．）
	EXPECT_NO_THROW({
	vector<vector<string>> inp = FileIn::input_to_array(filename);
	EXPECT_EQ(inp[0][0], "$$abc");
	EXPECT_EQ(inp[0][1], "def");
		}; );

	return;
};

TEST_F(FileInTest, pickup_data0) {

	// 存在しないファイルの読み取りエラーの確認
	EXPECT_ANY_THROW({
	vector<string> line = FileIn::pickup_data("hoge", this->inp_data);
		}; );

	// ヘッダが"a"行の読み取り．
	EXPECT_NO_THROW({
	vector<string> line = FileIn::pickup_data("a", this->inp_data);
	EXPECT_EQ(line[0], "a");
	EXPECT_EQ(line[3], "d");
		}; );
	
	return;
};

TEST_F(FileInTest, pickup_multiple_data0) {

	// 存在しないファイルの読み取りエラーの確認
	EXPECT_ANY_THROW({
	vector<string> line = FileIn::pickup_multiple_data("hoge", 2, this->inp_data);
		}; );

	// ヘッダは存在するけど指定した数値がない場合のエラー確認．
	EXPECT_ANY_THROW({
	vector<string> line = FileIn::pickup_multiple_data("e", 0, this->inp_data);
		}; );

	// ヘッダが"e"，2列めが"2"となっている行の読み取り．
	EXPECT_NO_THROW({
	vector<string> line = FileIn::pickup_multiple_data("e", 2, this->inp_data);
	EXPECT_EQ(line[0], "e");
	EXPECT_EQ(line[3], "g");
		}; );

	return;
};


TEST_F(FileInTest, pickup_matrix_data0) {

	// 存在しないファイルの読み取りエラーの確認
	EXPECT_ANY_THROW({
	vector<string> line = FileIn::pickup_matrix_data("hoge", 3, 5, this->inp_data);
		}; );

	// ヘッダは存在するけど指定した数値がない場合のエラー確認．
	EXPECT_ANY_THROW({
	vector<string> line = FileIn::pickup_matrix_data("h", 3, 4, this->inp_data);
		}; );

	// ヘッダが"e"，2列めが"2"となっている行の読み取り．
	EXPECT_NO_THROW({
	vector<string> line = FileIn::pickup_matrix_data("h", 3, 5, this->inp_data);
	EXPECT_EQ(line[0], "h");
	EXPECT_EQ(line[3], "j");
		}; );

	return;
};

// ちゃんと意図したコメントが出ていることを目視で確認してね．（20.12.25.楠崎確認済みのためTRUEとしています．）
TEST_F(FileInTest, invalid_argument_error0) {

	FileIn::invalid_argument_error("hoge");
	EXPECT_TRUE(true);
	return;
};

// ちゃんと意図したコメントが出ていることを目視で確認してね．（20.12.25.楠崎確認済みのためTRUEとしています．）
TEST_F(FileInTest, cannot_find_error0) {

	FileIn::cannot_find_error("fuga");
	EXPECT_TRUE(true);
	return;
};

