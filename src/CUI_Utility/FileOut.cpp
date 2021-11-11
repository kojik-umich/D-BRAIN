#include "FileOut.h"


// �t�@�C����V�K�쐬���C�w�b�_�[���������ރ��\�b�h�D
void FileOut::write_header(
	string name,	// [in]	�t�@�C����
	string head		// [in]	�������ރw�b�_�[������
) {
	ofstream file;
	file.open(name.c_str(), ios_base::trunc);
	file << head;

	return;
}

// �쐬����Ă���t�@�C���ɒǉ��Œl��1�s�ɏ������ރ��\�b�h�D
void FileOut::write_vector(
	string name,			// [in]	�t�@�C����
	const VectorXd&vec,		// [in]	�������ސ��l�z��
	int prec				// [in]	���l�̗L������
) {
	IOFormat CSVFormat(prec, 0, ", ", ", ", "", "", "", "\n");
	ofstream file;
	file.open(name.c_str(), ios_base::app);
	file << vec.format(CSVFormat);

	return;
}
