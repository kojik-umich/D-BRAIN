#pragma once

#include <iostream>
#include <fstream>
#include <Eigen\Dense>

using std::string;
using std::ofstream;
using std::ios_base;

using Eigen::IOFormat;
using Eigen::VectorXd;

namespace FileOut {
	void write_header(string name, string head);
	void write_vector(string name, const VectorXd & vector, int prec);
}
