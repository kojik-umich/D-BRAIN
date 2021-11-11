#include "FileIn.h"


// csv�t�@�C���̊e�s���x�N�g���ɓǂݍ��ނ��߂̃��\�b�h�Dhttps://www.sejuku.net/blog/48660 ����R�s�y
vector<string> FileIn::split
(					// out:	���͕��������؂蕶����(�܂ސ�:n)�ɂ����n+1�ɕ��������e������D
	string& input,	// in:	���͕�����D
	char delimiter	// in:	��؂蕶���D
) {
	istringstream stream(input);
	string field;
	vector<string> result;
	while (getline(stream, field, delimiter)) {
		result.push_back(field);
	}
	return result;
}

// csv�t�@�C���̂���$$�Ŏn�܂�s��string�^�z��Ɋi�[���郁�\�b�h
vector<vector<string>>FileIn::input_to_array
(								// out:	���̓t�@�C���S�̂�2�����z��D�i��v�f:�e�s�D���v�f:�e��j
	const char filename[]		// in:	���̓t�@�C�����D
) {
	//�t�@�C�����J���Ȃ�������G���[
	ifstream file(filename);
	if (!file) {
		cout << "Error! " << filename << " can not be opened" << endl;
		throw std::runtime_error("Inputfile Error");
	}

	// �e�s�ɑ΂��āC(0)$$�Ŏn�܂��2��ڂ���ǂݎ��(1)����ȊO�͓ǂݔ�΂��C�̃��[����(0)�ɊY���������̂����e�z��Ɋi�[�D
	string line;
	vector<vector<string>> inp_data;
	string dolmark = "$$";
	char delimiter = ',';
	while (getline(file, line)) {
		vector<string> strvec = FileIn::split(line, delimiter);
		size_t pos;

		//�X�y�[�X��S����������(���ӁF�ȉ��̃��[�v�����{����ƑS�p�����������������N����)
		while ((pos = (strvec.at(0)).find_first_of("�@ \t")) != string::npos) {
			strvec.at(0).erase(pos, 1);
		}

		//"$$"�Ŏn�܂�s�ł���Γǂݍ���
		if (strvec.at(0).size() >= dolmark.size() && equal(begin(dolmark), end(dolmark), begin(strvec.at(0)))) {
			inp_data.push_back(strvec);
		}
	}
	return inp_data;
}


// 1��ڂ�"param_name"���܂܂�Ă������������C���̍s�S�̂�string�z��Ɋi�[����D
vector<string> FileIn::pickup_data
(									// out:	�����ň�ԍŏ��Ƀq�b�g�����̍s���̓t�@�C��1�����z��D
	string param_name,				// in:	����������D
	const vector<vector<string>>&inp_data	// in:	���̓t�@�C���S�̂�2�����z��D�i��v�f:�e�s�D���v�f:�e��j
) {
	for (unsigned int i = 0; i < inp_data.size(); ++i)
		if (inp_data[i][0] == param_name)
			return inp_data[i];

	//$$~~��������Ȃ����G���[�𓊂���
	throw std::runtime_error(param_name + " is Not Found");
}

// 1��ڂ�"param_name"�C2��ڂ�elem_num���L�ڂ���Ă������������C���̍s�S�̂�string�z��Ɋi�[����D
vector<string> FileIn::pickup_multiple_data
(									// out:	�����ň�ԍŏ��Ƀq�b�g�����̍s���̓t�@�C��1�����z��D
	string param_name,				// in:	����������D
	int elem_num,					// in:	�v�f�ԍ��D
	const vector<vector<string>>&inp_data	// in:	���̓t�@�C���S�̂�2�����z��D�i��v�f:�e�s�D���v�f:�e��j
) {
	for (unsigned int i = 0; i < inp_data.size(); ++i)
		if (inp_data[i][0]  == param_name &&
			stoi(inp_data[i][1]) == elem_num)
			return inp_data[i];

	//$$~~��������Ȃ����G���[�𓊂���
	throw std::runtime_error("NotFound");
}

// 1��ڂ�"param_name"�C2��ڂ�"first_num"�C3��ڂ�"second_num"���L�ڂ���Ă������������C���̍s�S�̂�string�z��Ɋi�[����D
vector<string> FileIn::pickup_matrix_data
(									// out:	�����ň�ԍŏ��Ƀq�b�g�����̍s���̓t�@�C��1�����z��D
	string param_name,				// in:	����������D
	int first_num,					// in:	�v�f�ԍ��D
	int second_num,					// in:	�v�f�ԍ��D
	const vector<vector<string>>&inp_data	// in:	���̓t�@�C���S�̂�2�����z��D�i��v�f:�e�s�D���v�f:�e��j
) {
	for (unsigned int i = 0; i < inp_data.size(); ++i)
		if (inp_data[i][0]  == param_name &&
			stoi(inp_data[i][1]) == first_num &&
			stoi(inp_data[i][2]) == second_num)
			return inp_data[i];

	//$$~~��������Ȃ����G���[�𓊂���
	throw std::runtime_error("NotFound");
}

// ���̓t�@�C����$$hoge��̍s���ɕs�����������ۂɏo���G���[�D
bool FileIn::invalid_argument_error(const string & param_name) {
	cout <<  "�G���[�F�F" << param_name << "�̓��͂�����Ă��܂��D" << endl;
	cout <<  "�@�@�@�@�@input�t�@�C�����m�F���Ă��������D" << endl;
	return true;
}

// ���̓t�@�C����$$hoge�񂪖��������ۂɏo���G���[�D
bool FileIn::cannot_find_error(const string & param_name) {
	cout <<  "�G���[�F�F" << param_name << "�̓��͂�����܂���D" << endl;
	cout <<  "�@�@�@�@�@input�t�@�C�����m�F���Ă��������D" << endl;
	return true;
}



//while((pos = (strvec.at(0)).find_first_of(" ")) != string::npos){
			//if(strvec.at(0) == "$$"){
			//cout << strvec.at(0) << endl;
	//�ǂݎ��f�[�^�m�F�p
	//for(unsigned int i = 0; i < inp_data.size(); i++){
	//	for(unsigned int j =0; j < inp_data.at(i).size();j++){
	//		cout << inp_data.at(i).at(j) << ", ";
	//	}
	//	cout << endl;
	//}
