#include "pch.h"

class FileOutTest : public ::testing::Test {
public:
	string fname = "FileOut_Test.csv";
};

// 意図したテキストが出力されていることを目視にて確認．（20.12.25.楠崎）
TEST_F(FileOutTest, write_header_vector0) {

	FileOut::write_header(this->fname, "a,b,c,d,e,f,h\n");
	FileOut::write_vector(this->fname, VectorXd(7).setLinSpaced(2.0, 10.0), 10);

	EXPECT_TRUE(true);

	return;
};
