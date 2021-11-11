#pragma once

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;
using std::string;


namespace FileIn {

	vector<string> split(string& input, char delimiter);			// csv����擾�����������","�ŋ�؂��ăx�N�g���Ɋi�[���郁�\�b�h
	vector<vector<string>> input_to_array(const char filename[]);	// csv�t�@�C���S�̂�1��ڂ�string�^�z��Ɋi�[���郁�\�b�h
	vector<string> pickup_data(string param_name, const vector<vector<string>>&inp_data);
	vector<string> pickup_multiple_data(string param_name, int elem_num, const vector<vector<string>>&inp_data);
	vector<string> pickup_matrix_data(string param_name, int first_num, int second_num, const vector<vector<string>>&inp_data);
	bool invalid_argument_error(const string&param_name);
	bool cannot_find_error(const string&param_name);
};
