/*******************************************************************************
!								"B4P_FileIn.cpp"
!													2020/02/10	[Core-T]	���
!
!	�E4�_�ڐG�ʎ���̓���́E�É�͂��ׂĂ̓��͒l�̓ǂݍ��݂����̃N���X���ŒS��
!	�E���O�ցE�ێ���d�ʂ���O�ւ̍aR���S���a�E�ڐG�p�Ȃ�
!	�u�蓮���͂Ǝ����v�Z�̗������\�ȃp�����[�^�v�͎����v�Z�����̃N���X���Ŏ��{
!	�E�����v�Z�͕ʌɃ��\�b�h���쐬���C�v�Z�ɕK�v�Ȑ��l�͂��ׂĈ����ɐݒ�C
!	�@�Ԃ�l���v�Z���ʂƂ���D���̎����v�Z���s�����\�b�h�ł̓N���X�̃����o�ϐ��͎g�p�֎~
!	�E�w���\�L���ꂽ�����i1.0e+6�Ȃǁj��string�^����int�^�ɕϊ�����Ƃ���stoi()�ŕϊ��ł��Ȃ��̂ŁC
!	�@stod()���g�p���l�̌ܓ�����int�^�ɂ���D
!	�E�y�����K���z"read_~~~"�Ŏn�܂�֐��F�C���v�b�g�t�@�C������"$$~~~"�̓ǂݍ��݁C
!	�@�@�@�@�@�@�@"calc_~~"�Ŏn�܂�֐��F�e�퐔�l�̎����v�Z�D
!*******************************************************************************/

#include "B4P_FileIn.h"
using Numeric::Square;
using Numeric::Cube;

void B4P_FileIn::Nothing(void) {
}

// CSV�t�@�C����ǂݍ��݁C�e�����o�ϐ��Ɋi�[
bool B4P_FileIn::read_input_all(char fpath_d4bin_csv[]) {

	// �t�@�C����ǂݍ����string�^�̓񎟌��x�N�g���ɏ�������
	vector<vector<string>> inp_data;
	inp_data = FileIn::input_to_array(fpath_d4bin_csv);

	// �t�@�C���ǂݎ��G���[�̃t���b�O�itrue:�G���[����Cfalse:�G���[�����j
	bool has_error = false;

	// �􉽌`��̓ǂݎ��
	has_error |= this->read_RollingResistance(inp_data);
	has_error |= this->read_Coulomb(inp_data);
	has_error |= this->read_FilmThickness(inp_data);
	has_error |= this->read_Hysteresis(inp_data);


	// �]���̐��̎擾
	has_error |= this->read_BallNum(inp_data);
	int Z = this->ballnum;							// �����o�ϐ��̌Ăяo���͎��Ԃ�������̂ŁC���[�J���ϐ��ɏ����o��
	this->BL = new B4P_In::Ball[Z];
	for (int i = 0; i < Z; i++) {
		has_error |= this->read_Ball(inp_data, i);
	}
	
	has_error |= this->read_Inner(inp_data);
	has_error |= this->read_Outer(inp_data);
	this->ballpcd = (this->IR.Rod[0] + this->IR.Rod[1] + this->OR.Rod[0] + this->OR.Rod[1]) / 4;
	this->balldia = this->BL[0].dia;
	double cos_alpi0 = calc_cos_alp0(this->IR.R[0], this->balldia, this->IR.Rox[0]);
	double cos_alpi1 = calc_cos_alp0(this->IR.R[1], this->balldia, this->IR.Rox[1]);
	double cos_alpo0 = calc_cos_alp0(this->OR.R[0], this->balldia, this->OR.Rox[0]);
	double cos_alpo1 = calc_cos_alp0(this->OR.R[1], this->balldia, this->OR.Rox[1]);
	this->cos_alp0 = (cos_alpi0 + cos_alpi1 + cos_alpo0 + cos_alpo1) / 4;


	has_error |= this->read_SetCage(inp_data);

	switch (cage_type) {
	case snap_cage:
		has_error |= this->read_SnapCage(inp_data);
		break;

	default:
		cout << "�G���[�F���Ή��̕ێ��킪���͂���Ă��܂�" << endl;
		has_error = true;
		break;
	}

	has_error |= this->read_Oil(inp_data);
	has_error |= this->read_ContaBI(inp_data);
	has_error |= this->read_ContaBO(inp_data);
	has_error |= this->read_ContaBC(inp_data);

	has_error |= this->read_Gravity(inp_data);

	has_error |= this->read_LoadIn(inp_data);
	has_error |= this->read_InnnerBound(inp_data);
	has_error |= this->read_InnnerPosition(inp_data);
	has_error |= this->read_Dimension(inp_data);

	has_error |= this->read_Rotation(inp_data);
	has_error |= this->read_DynSet(inp_data);
	has_error |= this->read_FromContinuation(inp_data);
	has_error |= this->read_StopCalculation(inp_data);

	// �ϐ����̌v�Z
	this->DynSet.nX = 13 * ballnum + 39;
	this->SttSet[0].nX = 2 * ballnum + 5;
	this->SttSet[1].nX = 4 * ballnum;
	this->SttSet[2].nX = 6 * ballnum;

	// �f�t�H���g�̏o�̓t�@�C�������w��
	string fname_out = "d4b";
	has_error |= this->read_FileName(inp_data, fname_out);
	has_error |= this->read_OutputSlice(inp_data);
	return has_error;
}

//  BallNum�f�[�^�擾
bool B4P_FileIn::read_BallNum(const vector<vector<string>> &inp_data) {
	string param_name = "$$BallNum";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->ballnum		= stoi(param[1]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}


//  SetCage�f�[�^�擾
bool B4P_FileIn::read_SetCage(const vector<vector<string>> &inp_data) {
	string param_name = "$$SetCage";
	try {
		vector<string> param	= FileIn::pickup_data(param_name, inp_data);
		this->cage_type			= static_cast<CageType>(stoi(param[1]));
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}



//  Ball�f�[�^�擾
bool B4P_FileIn::read_Ball(const vector<vector<string>> &inp_data, int i) {
	string param_name = "$$Ball";
	try {
		for (int i = 0; i < this->ballnum; i++) {
			vector<string> param	= FileIn::pickup_multiple_data(param_name, i, inp_data);
			this->BL[i].den		= stod(param[2]);
			this->BL[i].E		= stod(param[3]);
			this->BL[i].por		= stod(param[4]);
			this->BL[i].rms		= stod(param[5]);
			this->BL[i].dia		= stod(param[6]);
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

//  Inner�f�[�^�擾
bool B4P_FileIn::read_Inner(const vector<vector<string>> &inp_data) {
	string param_name = "$$Inner";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->IR.den	= stod(param[1]);
		this->IR.E		= stod(param[2]);
		this->IR.por	= stod(param[3]);
		this->IR.rms	= stod(param[4]);
		this->IR.R[0]	= stod(param[5]);
		this->IR.R[1]	= stod(param[6]);
		this->IR.Rox[0]	= stod(param[7]);
		this->IR.Rox[1]	= stod(param[8]);
		this->IR.Rod[0] = stod(param[9]);
		this->IR.Rod[1] = stod(param[10]);
		this->IR.hedge[0]	= stod(param[11]);
		this->IR.hedge[1]	= stod(param[12]);

		this->IR.m = stod(param[13]);
		this->IR.Ix = stod(param[14]);
		this->IR.Iyz = stod(param[15]);

	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;

}
//  Outer�f�[�^�擾
bool B4P_FileIn::read_Outer(const vector<vector<string>> &inp_data) {
	string param_name = "$$Outer";
	vector<string> param = FileIn::pickup_data(param_name, inp_data);
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->OR.den		= stod(param[1]);
		this->OR.E			= stod(param[2]);
		this->OR.por		= stod(param[3]);
		this->OR.rms		= stod(param[4]);
		this->OR.R[0]		= stod(param[5]);
		this->OR.R[1]		= stod(param[6]);
		this->OR.Rox[0]		= stod(param[7]);
		this->OR.Rox[1]		= stod(param[8]);
		this->OR.Rod[0] = stod(param[9]);
		this->OR.Rod[1] = stod(param[10]);
		this->OR.hedge[0]	= stod(param[11]);
		this->OR.hedge[1]	= stod(param[12]);
		this->OR.m = stod(param[13]);
		this->OR.Ix = stod(param[14]);
		this->OR.Iyz = stod(param[15]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// ���^�ێ���e�p�����[�^�̎擾�E�ݒ�
bool B4P_FileIn::read_SnapCage(const vector<vector<string>> &inp_data) {

	string param_name = "$$SnapCage";

	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->Cage.den	= stod(param[1]);
		this->Cage.E	= stod(param[2]);
		this->Cage.por	= stod(param[3]);
		this->Cage.rms	= stod(param[4]);
		this->Snap.dout	= stod(param[5]);
		this->Snap.din	= stod(param[6]);
		double clearance = stod(param[7]);		// �|�P�b�g����
		this->Snap.R    = 0.5 * (this->balldia + clearance);
		this->Snap.h	= stod(param[8]);
		this->Snap.jc	= stod(param[9]);

		this->Cage.rmg0[0] = stod(param[10]);
		this->Cage.rmg0[1] = stod(param[11]);
		this->Cage.rmg0[2] = stod(param[12]);
		this->Snap.ropen	= stod(param[13]);
		this->Cage.m = stod(param[14]);
		this->Cage.Ix = stod(param[15]);
		this->Cage.Iyz = stod(param[16]);

		// �|�P�b�g���̃p�����[�^��`
		this->Snap.Kface		= stod(param[17]);
		this->Snap.Kedgein		= stod(param[18]);
		this->Snap.Kedgeout		= stod(param[19]);
		this->Snap.Kcornerin	= stod(param[20]);
		this->Snap.Kcornerout	= stod(param[21]);
		this->Snap.Kopen		= stod(param[22]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

//  Oil�f�[�^�擾
bool B4P_FileIn::read_Oil(const vector<vector<string>> &inp_data) {

	string param_name = "$$Oil";

	try {
		for (int i = 0; i < this->ballnum; i++) {
			vector<string> param = FileIn::pickup_data(param_name, inp_data);
			this->LB.alpha0	= Unit::mm2kgf2Painv(stod(param[1]));
			this->LB.beta0	= stod(param[2]);
			this->LB.k0		= stod(param[3]);
			this->LB.eta0	= stod(param[4]);
			this->LB.lm0	= stod(param[5]);
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

//  ContaBI�f�[�^�擾
bool B4P_FileIn::read_ContaBI(const vector<vector<string>> &inp_data) {

	string param_name = "$$ContaBI";

	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->BIP.mu	= stod(param[1]);
		this->BIP.dzeta	= stod(param[2]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

//  ContaBO�f�[�^�擾
bool B4P_FileIn::read_ContaBO(const vector<vector<string>> &inp_data) {

	string param_name = "$$ContaBO";

	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->BOP.mu	= stod(param[1]);
		this->BOP.dzeta	= stod(param[2]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

//  ContaBC�f�[�^�擾
bool B4P_FileIn::read_ContaBC(const vector<vector<string>> &inp_data) {
	string param_name = "$$ContaBC";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->BCP.mu	= stod(param[1]);
		this->BCP.dzeta	= stod(param[2]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

//  ContaCI�f�[�^�擾
bool B4P_FileIn::read_ContaCI(const vector<vector<string>> &inp_data) {

	string param_name = "$$ContaCI";

	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->CIP.mu	= stod(param[1]);
		this->CIP.dzeta	= stod(param[2]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

//  ContaCO�f�[�^�擾
bool B4P_FileIn::read_ContaCO(const vector<vector<string>> &inp_data) {

	string param_name = "$$ContaCO";

	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->COP.mu	= stod(param[1]);
		this->COP.dzeta	= stod(param[2]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

//  Gravity�f�[�^�擾
bool B4P_FileIn::read_Gravity(const vector<vector<string>> &inp_data) {
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

//  LoadIn�f�[�^�擾
bool B4P_FileIn::read_LoadIn(const vector<vector<string>> &inp_data) {

	string param_name = "$$LoadIn";

	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->LoadIn[0]			= stod(param[1]);
		this->LoadIn[1]			= stod(param[2]);
		this->LoadIn[2]			= stod(param[3]);
		this->LoadIn[3]			= 0;
		this->LoadIn[4]			= stod(param[4]);
		this->LoadIn[5]			= stod(param[5]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

//  Rotation�f�[�^�擾
bool B4P_FileIn::read_Rotation(const vector<vector<string>> &inp_data) {

	string param_name = "$$Rotation";

	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);

		this->omegair			= Unit::rpm2rads(stod(param[1]));
		this->omegaor			= Unit::rpm2rads(stod(param[2]));
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

//  Dynamic�f�[�^�擾
bool B4P_FileIn::read_DynSet(const vector<vector<string>> &inp_data) {

	string param_name = "$$DynSet";

	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->DynSet.calctime	= stod(param[1]);
		this->DynSet.sampling	= Numeric::Roundoff(stod(param[2]));
		this->DynSet.h			= stod(param[3]);
		this->DynSet.hmin		= stod(param[4]);
		this->DynSet.ep		= stod(param[5]);
		this->DynSet.tr		= stod(param[6]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

//  Static�f�[�^�擾(�g���ĂȂ�)
bool B4P_FileIn::read_SttSet(const vector<vector<string>> &inp_data, int i) {
	string param_name = "$$SttSet" + to_string(i);
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->SttSet[i].iter1 = stoi(param[1]);
		this->SttSet[i].iter2 = stoi(param[2]);
		this->SttSet[i].rs = stod(param[3]);
		this->SttSet[i].jac_eps = stod(param[4]);
		for (int j = 0; j < 6; j++)
			this->SttSet[i].eps[j] = stod(param[j + 5]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

//  Dynamic�ψʍS�������f�[�^�擾
bool B4P_FileIn::read_InnnerPosition(const vector<vector<string>> &inp_data) {
	string param_name = "$$InnnerPosition";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->DynSet.x[0]	= Unit::mm2m(stod(param[1]));
		this->DynSet.x[1]	= Unit::mm2m(stod(param[2]));
		this->DynSet.x[2]	= Unit::mm2m(stod(param[3]));
		this->DynSet.ax[0] = stod(param[4]);
		this->DynSet.ax[1] = stod(param[5]);
		this->DynSet.ax[2] = stod(param[6]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}



//  Dynamic�S�������f�[�^�擾
bool B4P_FileIn::read_InnnerBound(const vector<vector<string>> &inp_data) {

	string param_name = "$$InnnerBound";
	bool has_error = false;
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->bound.v_const[0]	= (stoi(param[1]) > 0);
		this->bound.v_const[1]	= (stoi(param[2]) > 0);
		this->bound.v_const[2]	= (stoi(param[3]) > 0);
		this->bound.w_const[0]	= (stoi(param[4]) > 0);
		this->bound.w_const[1]	= (stoi(param[5]) > 0);
		this->bound.w_const[2]	= (stoi(param[6]) > 0);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return has_error;
}

//  FileName�f�[�^�擾�i���͒l�Ɍ�肪����ꍇ�C�f�t�H���g�̃t�@�C������ݒ�j
bool B4P_FileIn::read_FileName(const vector<vector<string>> &inp_data, string fname_out) {

	string param_name = "$$FileName";

	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->FN.Temp			= param[1];
		this->FN.Ball			= param[2];
		this->FN.Outer			= param[3];
		this->FN.Inner			= param[4];
		this->FN.Cage			= param[5];
		this->FN.BallOuterPair	= param[6];
		this->FN.BallInnerPair	= param[7];
		this->FN.BallCagePair	= param[8];
	}
	catch (invalid_argument) {
		cout <<  "���ӁF�F" << param_name << "�̓��͂�����Ă��܂��D" << endl;
		cout <<  "�@�@�@�@�@�f�t�H���g�̃t�@�C������ݒ肵�܂��D" << endl;
		this->FN.Temp			= "00_Temp";
		this->FN.Ball			= "01_Ball";
		this->FN.Outer			= "02_Outer";
		this->FN.Inner			= "03_Inner";
		this->FN.Cage			= "04_Cage";
		this->FN.BallOuterPair	= "05_BallOuterPair";
		this->FN.BallInnerPair	= "06_BallInnerPair";
		this->FN.BallCagePair	= "07_BallCagePair";
	}
	catch (exception) {
		cout <<  "���ӁF�F" << param_name << "�̓��͂�����܂���D" << endl;
		cout <<  "�@�@�@�@�@�f�t�H���g�̃t�@�C������ݒ肵�܂��D" << endl;
		this->FN.Temp			= "_00_Temp";
		this->FN.Ball			= "_01_Ball";
		this->FN.Outer			= "_02_Outer";
		this->FN.Inner			= "_03_Inner";
		this->FN.Cage			= "_04_Cage";
		this->FN.BallOuterPair	= "_05_BallOuterPair";
		this->FN.BallInnerPair	= "_06_BallInnerPair";
		this->FN.BallCagePair	= "_07_BallCagePair";
	}
	return false;
}

//  �]���薀�C�ݒ�����f�[�^�擾
bool B4P_FileIn::read_RollingResistance(const vector<vector<string>> &inp_data) {

	string param_name = "$$RollingResistance";

	try {
		vector<string> param	= FileIn::pickup_data(param_name, inp_data);
		this->TB.rollingresistance = static_cast<Tribology::RollingResistance>(stoi(param[1]));
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// �����ڐG���薀�C�ݒ�����f�[�^�擾
bool B4P_FileIn::read_Coulomb(const vector<vector<string>>&inp_data) {
	string param_name = "$$Coulomb";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->TB.coulomb = static_cast<Tribology::Coulomb>(stoi(param[1]));
		this->TB.coulomb_slope = stod(param[2]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// �q�X�e���V�X�����f�[�^�擾
bool B4P_FileIn::read_Hysteresis(const vector<vector<string>>&inp_data) {
	string param_name = "$$Hysteresis";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->TB.hysteresis = static_cast<Tribology::Hysteresis>(stoi(param[1]));
		this->TB.hysteresis_factor = stod(param[2]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// ���������v�Z���f�[�^�擾
bool B4P_FileIn::read_FilmThickness(const vector<vector<string>>&inp_data) {
	string param_name = "$$FilmThickness";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->TB.filmThickness = static_cast<Tribology::FilmThickness>(stoi(param[1]));
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

bool B4P_FileIn::read_Dimension(const vector<vector<string>> &inp_data) {

	string param_name = "$$Dimension";

	try {
		vector<string> param	= FileIn::pickup_data(param_name, inp_data);
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

// �X���C�X�o�͐ݒ�
bool B4P_FileIn::read_OutputSlice(const vector<vector<string>> &inp_data) {

	string param_name = "$$OutputSlice";

	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		int n = stoi(param[1]);
		this->OutSlice.n = n;
		this->OutSlice.list = new(int[n]);
		for (int i = 0; i < n; i++) {
			this->OutSlice.list[i] = stoi(param[i + 2]);
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

// ����͑����v�Z
bool B4P_FileIn::read_FromContinuation(const vector<vector<string>> &inp_data) {
	string param_name = "$$FromContinuation";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->DynSet.fromcontinuation = (stoi(param[1]) == 1);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}

// ����͌v�Z�ł��؂�
bool B4P_FileIn::read_StopCalculation(const vector<vector<string>> &inp_data) {
	string param_name = "$$StopCalculation";
	try {
		vector<string> param = FileIn::pickup_data(param_name, inp_data);
		this->DynSet.stopcalculation = (stoi(param[1]) == 1);
		this->DynSet.dTierr = stod(param[2]);
		this->DynSet.stp = stoi(param[3]);
	}
	catch (invalid_argument) {
		return FileIn::invalid_argument_error(param_name);
	}
	catch (exception) {
		return FileIn::cannot_find_error(param_name);
	}
	return false;
}
// �ڐG�p�̌v�Z�i�ʓx�@�ŏo�́j
double B4P_FileIn::calc_cos_alp0(
	double R,		//  �aR[m]
	double D,		//  �ʌa[m] 
	double RO_x		//  �aR���Sx���W[m]
) {
	return  sqrt(Square(R - D / 2) - Square(RO_x)) / (R - D / 2);
}



















////  AnalysisMode�f�[�^�擾
//bool B4P_FileIn::read_AnalysisMode(const vector<vector<string>> &inp_data) {
//	string param_name = "$$AnalysisMode";
//	try {
//		vector<string> param = FileIn::pickup_data(param_name, inp_data);
//		this->AnaMode = static_cast<AnalysisMode>(stoi(param[1]));
//	}
//	catch (invalid_argument) {
//		return FileIn::invalid_argument_error(param_name);
//	}
//	catch (exception) {
//		return FileIn::cannot_find_error(param_name);
//	}
//	return false;
//}












// for(int i = 0; i < this->ballnum; i++){
// 	this->Ball[i].DInit("$$Ball", 6, true);
// 	this->Ball[i].param		=  this->multi_dparam_read(this->Ball[i], inp_data, (i + 1), this->ballnum);
// }
// this->Inner.DInit("$$Inner", 13, true);
// this->Outer.DInit("$$Outer", 13, true);
// �d�ʎ����v�Z
// if(this->Outer.state[12] == auto_calc)
// 	this->Outer.param[12] = calc_ringmass(this->ballpcd, Outer.param[0]);
// if(this->Inner.state[12] == auto_calc)
// 	this->Inner.param[12] = calc_ringmass(this->ballpcd, Inner.param[0]);
// cout << "�O�֏d�ʎ����v�Z" << Outer.param[12] << endl;
// cout << "���֏d�ʎ����v�Z" << Inner.param[12] << endl;
//  has_error |= this->read_ContaCI(inp_data);
//  has_error |= this->read_ContaCO(inp_data);
//  has_error |= this->read_DisplIn(inp_data);
//  has_error |= this->read_DisplOut(inp_data);
//  has_error |= this->read_DisplCage(inp_data);
//  has_error |= this->read_InitBall(inp_data);
// dbb�ɏ����銥�^�ێ��펿�ʂ̌v�Z
// double B4P_FileIn::calc_Snap_mass(double Snap_den_, double Snap_rout_, double Snap_rin_, double Cage_h_, double ballpcd_, double balldia_, double clr_, double ballnum_){
// 	double rp3 = 0.5*ballpcd_+0.5*balldia_+0.25*clr_;// �O�֔��a
// 	double mc;
// 	/* snap cage 1.4465 is modification factor */
// 	mc=1.4465*Snap_den_*Numeric::pi*((Snap_rout_*Snap_rout_-Snap_rin_*Snap_rin_)*Cage_h_
//     	-(double)ballnum_*(Snap_rout_-Snap_rin_)*rp3*rp3
// 		*(0.5+0.5*sin(2.0*asin((Cage_h_)/rp3))+asin((Cage_h_)/rp3)));
// 	return mc;
// }
		//  �ē������a	= stod(param[13])
		//  �ē�����		= stod(param[14])
		//  �ē��ʑ��p	= stod(param[15])

			// this->autocalc_Inner_Iy = true;
			// this->autocalc_Inner_Iy = false;


// //  SigmLoad�f�[�^�擾
// bool B4P_FileIn::read_SigmLoad(const vector<vector<string>> &inp_data) {
// 
// 	string param_name = "$$SigmLoad";
// 
// 	try {
// 		vector<string> param = FileIn::pickup_data(param_name, inp_data);
// 		this->Time_Fstart		= stod(param[1]);
// 		this->Time_Fend			= stod(param[2]);
// 		this->Time_Rstart		= stod(param[3]);
// 		this->Time_Rend			= stod(param[4]);
// 		this->Sigm_srate		= stod(param[5]);
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }
// 
// //  LinearLoad�f�[�^�擾
// bool B4P_FileIn::read_LinearLoad(const vector<vector<string>> &inp_data) {
// 
// 	string param_name = "$$LinearLoad";
// 
// 	try {
// 		vector<string> param = FileIn::pickup_data(param_name, inp_data);
// 		this->Time_Fstart		= stod(param[1]);
// 		this->Time_Fend			= stod(param[2]);
// 		this->Time_Rstart		= stod(param[3]);
// 		this->Time_Rend			= stod(param[4]);
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }
// 
// //  LogLoad�f�[�^�擾
// bool B4P_FileIn::read_LogLoad(const vector<vector<string>> &inp_data) {
// 
// 	string param_name = "$$LogLoad";
// 
// 	try {
// 		vector<string> param = FileIn::pickup_data(param_name, inp_data);
// 		this->Time_Fstart		= stod(param[1]);
// 		this->Time_Fend			= stod(param[2]);
// 		this->Time_Rstart		= stod(param[3]);
// 		this->Time_Rend			= stod(param[4]);
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }
// //  SetLoad�f�[�^�擾
// bool B4P_FileIn::read_SetLoad(const vector<vector<string>> &inp_data) {
// 	string param_name = "$$SetLoad";
// 	try {
// 		vector<string> param	= FileIn::pickup_data(param_name, inp_data);
// 		this->load_change		= static_cast<IncreaseFuncType>(stoi(param[1]));;
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }
	// switch (this->load_change) {
	// case expotential_increase:
	// 	has_error |= this->read_LogLoad(inp_data);
	// 	break;

	// case constant_value:
	// 	break;

	// case sigmoid_increase:
	// 	has_error |= this->read_SigmLoad(inp_data);
	// 	break;

	// case linear_increase:
	// 	has_error |= this->read_LinearLoad(inp_data);
	// 	break;

	// case user_define:
	// default:
	// 	cout << "���ӁF$$SetLoad���K�؂ɓ��͂���Ă��Ȃ����߁C"		<< endl
	// 		<< "�@�@�@�׏d�͒萔�l�Ƃ��ēǂݍ��݂܂�"				<< endl;
	// 	break;
	// }
// 
// //  SetStr�f�[�^�擾
// bool B4P_FileIn::read_SetStr(const vector<vector<string>> &inp_data) {
// 	string param_name = "$$SetStr";
// 	try {
// 		vector<string> param	= FileIn::pickup_data(param_name, inp_data);
// 		this->setstr			=  static_cast<StrType>(stoi(param[1]));
// 
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }

// //  DisplIn�f�[�^�擾
// bool B4P_FileIn::read_DisplIn(const vector<vector<string>> &inp_data) {
// 
// 	string param_name = "$$DisplIn";
// 
// 	try {
// 		vector<string> param = FileIn::pickup_data(param_name, inp_data);
// 		this->DisplIn[0]			= (stoi(param[1]) == 1);
// 		this->DisplIn[1]			= (stoi(param[2]) == 1);
// 		this->DisplIn[2]			= (stoi(param[3]) == 1);
// 		this->DisplIn[3]			= (stoi(param[4]) == 1);
// 		this->DisplIn[4]			= (stoi(param[5]) == 1);
// 		this->DisplIn[5]			= (stoi(param[6]) == 1);
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }
// //  DisplOut�f�[�^�擾
// bool B4P_FileIn::read_DisplOut(const vector<vector<string>> &inp_data) {
// 
// 	string param_name = "$$DisplOut";
// 
// 	try {
// 		vector<string> param = FileIn::pickup_data(param_name, inp_data);
// 		this->DisplOut[0]			= (stoi(param[1]) == 1);
// 		this->DisplOut[1]			= (stoi(param[2]) == 1);
// 		this->DisplOut[2]			= (stoi(param[3]) == 1);
// 		this->DisplOut[3]			= (stoi(param[4]) == 1);
// 		this->DisplOut[4]			= (stoi(param[5]) == 1);
// 		this->DisplOut[5]			= (stoi(param[6]) == 1);
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }
// //  DisplCage�f�[�^�擾
// bool B4P_FileIn::read_DisplCage(const vector<vector<string>> &inp_data) {
// 
// 	string param_name = "$$DisplCage";
// 
// 	try {
// 		vector<string> param = FileIn::pickup_data(param_name, inp_data);
// 		this->DisplCage[0] = (stoi(param[1]) != 0);
// 		this->DisplCage[1] = (stoi(param[2]) != 0);
// 		this->DisplCage[2] = (stoi(param[3]) != 0);
// 		this->DisplCage[3] = (stoi(param[4]) != 0);
// 		this->DisplCage[4] = (stoi(param[5]) != 0);
// 		this->DisplCage[5] = (stoi(param[6]) != 0);
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }
// 
// //  LoadOut�f�[�^�擾
// bool B4P_FileIn::read_LoadOut(const vector<vector<string>> &inp_data) {
// 
// 	string param_name = "$$LoadOut";
// 
// 	try {
// 		vector<string> param = FileIn::pickup_data(param_name, inp_data);
// 		this->LoadOut			= VectorXd::Zero(6);
// 		this->LoadOut[0]		= stod(param[1]);
// 		this->LoadOut[1]		= stod(param[2]);
// 		this->LoadOut[2]		= stod(param[3]);
// 		this->LoadOut[3]		= 0;
// 		this->LoadOut[4]		= stod(param[4]);
// 		this->LoadOut[5]		= stod(param[5]);
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }
// 
// //  LoadCage�f�[�^�擾
// bool B4P_FileIn::read_LoadCage(const vector<vector<string>> &inp_data) {
// 
// 	string param_name = "$$LoadCage";
// 
// 	try {
// 		vector<string> param = FileIn::pickup_data(param_name, inp_data);
// 		this->LoadCage			= VectorXd::Zero(6);
// 		this->LoadCage[0]		= stod(param[1]);
// 		this->LoadCage[1]		= stod(param[2]);
// 		this->LoadCage[2]		= stod(param[3]);
// 		this->LoadCage[3]		= 0;
// 		this->LoadCage[4]		= stod(param[4]);
// 		this->LoadCage[5]		= stod(param[5]);
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }
// //  InitIn�f�[�^�擾
// bool B4P_FileIn::read_InitIn(const vector<vector<string>> &inp_data) {
// 
// 	string param_name = "$$InitIn";
// 
// 	try {
// 		vector<string> param = FileIn::pickup_data(param_name, inp_data);
// 		this->InitIn		= VectorXd::Zero(6);
// 		this->InitIn[0]		= stod(param[1]);
// 		this->InitIn[1]		= stod(param[2]);
// 		this->InitIn[2]		= stod(param[3]);
// 		this->InitIn[3]		= stod(param[4]);
// 		this->InitIn[4]		= stod(param[5]);
// 		this->InitIn[5]		= stod(param[6]);
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }
// //  InitCage�f�[�^�擾
// bool B4P_FileIn::read_InitOut(const vector<vector<string>> &inp_data) {
// 
// 	string param_name = "$$InitOut";
// 
// 	try {
// 		vector<string> param = FileIn::pickup_data(param_name, inp_data);
// 		this->InitOut			= VectorXd::Zero(6);
// 		this->InitOut[0]		= stod(param[1]);
// 		this->InitOut[1]		= stod(param[2]);
// 		this->InitOut[2]		= stod(param[3]);
// 		this->InitOut[3]		= stod(param[4]);
// 		this->InitOut[4]		= stod(param[5]);
// 		this->InitOut[5]		= stod(param[6]);
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }
// 
// //  InitCage�f�[�^�擾
// bool B4P_FileIn::read_InitCage(const vector<vector<string>> &inp_data) {
// 
// 	string param_name = "$$InitCage";
// 
// 	try {
// 		vector<string> param = FileIn::pickup_data(param_name, inp_data);
// 		this->InitCage			= VectorXd::Zero(6);
// 		this->InitCage[0]		= stod(param[1]);
// 		this->InitCage[1]		= stod(param[2]);
// 		this->InitCage[2]		= stod(param[3]);
// 		this->InitCage[3]		= stod(param[4]);
// 		this->InitCage[4]		= stod(param[5]);
// 		this->InitCage[5]		= stod(param[6]);
// 	}
// 	catch (invalid_argument) {
// 		return FileIn::invalid_argument_error(param_name);
// 	}
// 	catch (exception) {
// 		return FileIn::cannot_find_error(param_name);
// 	}
// 	return false;
// }

