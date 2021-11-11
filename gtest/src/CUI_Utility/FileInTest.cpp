#include "pch.h"

class FileInTest : public ::testing::Test {
public:
	vector<vector<string>> inp_data = {
		{"a", "b", "c", "d"} , {"e", "2", "f", "g"}, {"h", "3", "5", "j"} };
};

TEST_F(FileInTest, input_to_array0) {

	// ���݂��Ȃ��t�@�C���̓ǂݎ��G���[�̊m�F
	EXPECT_ANY_THROW({
	vector<vector<string>> inp = FileIn::input_to_array("hoge.csv");
		}; );

	// out�ɒu���Ă���t�H���_���Q�ƁD
	char filename[] = "FileIn_test.csv";

	// �e�L�X�g�t�@�C���̓ǂݍ��݁D�i�G���[���o�����proj�̃v���p�e�B���ǂ����DOpenCppCoverage���Y�ꂸ�ɁD�j
	EXPECT_NO_THROW({
	vector<vector<string>> inp = FileIn::input_to_array(filename);
	EXPECT_EQ(inp[0][0], "$$abc");
	EXPECT_EQ(inp[0][1], "def");
		}; );

	return;
};

TEST_F(FileInTest, pickup_data0) {

	// ���݂��Ȃ��t�@�C���̓ǂݎ��G���[�̊m�F
	EXPECT_ANY_THROW({
	vector<string> line = FileIn::pickup_data("hoge", this->inp_data);
		}; );

	// �w�b�_��"a"�s�̓ǂݎ��D
	EXPECT_NO_THROW({
	vector<string> line = FileIn::pickup_data("a", this->inp_data);
	EXPECT_EQ(line[0], "a");
	EXPECT_EQ(line[3], "d");
		}; );
	
	return;
};

TEST_F(FileInTest, pickup_multiple_data0) {

	// ���݂��Ȃ��t�@�C���̓ǂݎ��G���[�̊m�F
	EXPECT_ANY_THROW({
	vector<string> line = FileIn::pickup_multiple_data("hoge", 2, this->inp_data);
		}; );

	// �w�b�_�͑��݂��邯�ǎw�肵�����l���Ȃ��ꍇ�̃G���[�m�F�D
	EXPECT_ANY_THROW({
	vector<string> line = FileIn::pickup_multiple_data("e", 0, this->inp_data);
		}; );

	// �w�b�_��"e"�C2��߂�"2"�ƂȂ��Ă���s�̓ǂݎ��D
	EXPECT_NO_THROW({
	vector<string> line = FileIn::pickup_multiple_data("e", 2, this->inp_data);
	EXPECT_EQ(line[0], "e");
	EXPECT_EQ(line[3], "g");
		}; );

	return;
};


TEST_F(FileInTest, pickup_matrix_data0) {

	// ���݂��Ȃ��t�@�C���̓ǂݎ��G���[�̊m�F
	EXPECT_ANY_THROW({
	vector<string> line = FileIn::pickup_matrix_data("hoge", 3, 5, this->inp_data);
		}; );

	// �w�b�_�͑��݂��邯�ǎw�肵�����l���Ȃ��ꍇ�̃G���[�m�F�D
	EXPECT_ANY_THROW({
	vector<string> line = FileIn::pickup_matrix_data("h", 3, 4, this->inp_data);
		}; );

	// �w�b�_��"e"�C2��߂�"2"�ƂȂ��Ă���s�̓ǂݎ��D
	EXPECT_NO_THROW({
	vector<string> line = FileIn::pickup_matrix_data("h", 3, 5, this->inp_data);
	EXPECT_EQ(line[0], "h");
	EXPECT_EQ(line[3], "j");
		}; );

	return;
};

// �����ƈӐ}�����R�����g���o�Ă��邱�Ƃ�ڎ��Ŋm�F���ĂˁD�i20.12.25.���m�F�ς݂̂���TRUE�Ƃ��Ă��܂��D�j
TEST_F(FileInTest, invalid_argument_error0) {

	FileIn::invalid_argument_error("hoge");
	EXPECT_TRUE(true);
	return;
};

// �����ƈӐ}�����R�����g���o�Ă��邱�Ƃ�ڎ��Ŋm�F���ĂˁD�i20.12.25.���m�F�ς݂̂���TRUE�Ƃ��Ă��܂��D�j
TEST_F(FileInTest, cannot_find_error0) {

	FileIn::cannot_find_error("fuga");
	EXPECT_TRUE(true);
	return;
};

