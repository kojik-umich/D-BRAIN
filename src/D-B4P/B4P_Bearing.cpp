/*******************************************************************************
!								"B4P_Bearing.cpp"
!													2018/12/18	[Core-T]	���
!
!	4�_�ڐG�ʎ��󂻂̂��̂̃I�u�W�F�N�g�D
!
!	���݂̏�Ԃł̑��݉׏d��O���[�o���ϐ����C�N�Ȓl�̕ێ��ȂǁC
!	���\�ȐӖ��𕉂��Ă����ςȐl�D�i�_�N���X�ɂȂ�Ȃ��悤�ɒ��ӁD�j
!	���M���ׂ��_�Ƃ���B4P_In�I�u�W�F�N�g���o�����Ēl������Ӗ��͂�����1�l�ŕ����Ă���D
!
!*******************************************************************************/

#include "B4P_Bearing.h"

// �v�Z�̓��͏�����ǂݍ��ށD�O�����K�v�ȃf�[�^�͓n���Ă��D
void B4P_Bearing::init(const B4P_In&FI, double*x, double*ax) {

	// �e���ނ̃W�I���g���ƍޗ��̐ݒ�D
	double Cage_rmg0[3];
	this->init_parts(FI, Cage_rmg0);

	// �����ʒu�C�������x��������D
	if (true)
		this->init_position(FI.balldia, FI.ballnum, FI.ballpcd, FI.cos_alp0, FI.omegair, FI.omegaor, Cage_rmg0);

	// ���̓t�@�C�����珉���ʒu�����肷��D(������)
	else
		double hoge;

	//	�ψʁE���x�̎��R�x�Ɨ́E���[�����g�̎��R�x�̕������m�ہD
	this->nX = 13 * FI.ballnum + 26;

	// ���ւ̍S��������ݒ�D
	this->IR.set_const(FI.bound.v_const, FI.bound.w_const);
	// ���ւ̏����l��ݒ�
	this->IR.x = Vector3d(x);
	this->IR.set_ax(Vector3d(ax));
	double w = FI.omegair;
	this->IR.w = w * this->IR.get_ax();

	return;
}

void B4P_Bearing::allocate_Pair(int Z) {

	this->BOP = new B4P_BallOuterRingPair[Z];
	this->BIP = new B4P_BallInnerRingPair[Z];
	this->BCP = new B4P_BallCagePair[Z];

	return;
}

// B4P_In����K�v�Ȓl��b4p�ɓ���郁�\�b�h�D�W�I���g���ƍޗ��D
void B4P_Bearing::init_parts(const B4P_In&FI, double*Cage_rmg0) {

	// �{�[�����������o���C�z��̓��I�m�ہD
	this->Z = FI.ballnum;
	this->BL = new Ball[Z];

	// �e���ނ̃W�I���g���̐ݒ�D
	this->OR.init(FI.OR, FI.ballpcd);
	this->IR.init(FI.IR, FI.ballpcd);
	switch (FI.cage_type) {
	case B4P_In::snap_cage:
		this->CG = new bal_SnapCage;
		for (int i = 0; i < 3; i++)
			Cage_rmg0[i] = FI.Cage.rmg0[i];
		break;
	default:
		this->CG = new bal_SnapCage;
		break;
	}
	this->CG->init(FI);

	bool v_const[3], w_const[3];
	for (int i = 0; i < 3; i++) {
		v_const[i] = false;
		w_const[i] = false;
	}

	for (int i = 0; i < this->Z; i++)
		this->BL[i].init(FI.BL[i].dia, FI.BL[i].E, FI.BL[i].por, FI.BL[i].den, FI.BL[i].rms, v_const, w_const);

	this->allocate_Pair(this->Z);

	// ball_X_pair�̑����R�Â���D
	for (int i = 0; i < this->Z; i++) {
		this->BIP[i].link(&this->BL[i], &this->IR);
		this->BOP[i].link(&this->BL[i], &this->OR);
		this->BCP[i].link(&this->BL[i], this->CG, i);
	}

	// ball_X_pair�̏��������s���D
	for (int i = 0; i < this->Z; i++) {
		this->BIP[i].init(FI);
		this->BOP[i].init(FI);
		this->BCP[i].init(FI);
	}

	// ���ɕK�v�ȏ���b4p�ɕێ�����D
	this->pcd = FI.ballpcd;

	Vector3d Fin = Vector3d(FI.LoadIn[0], FI.LoadIn[1], FI.LoadIn[2]);
	Vector3d Nin = Vector3d(FI.LoadIn[3], FI.LoadIn[4], FI.LoadIn[5]);

	// ����D
	this->F_load = Fin;
	this->T_load = Nin;

	// �e���ނɉ׏d�x�N�g����ݒ�D
	Rigid::g = Vector3d(FI.rigid.g);

	return;
}


// B4P_In���珉���ψʁE���x���Z�o���郁�\�b�h�D�v�Z�͒����ʒu�C���_���]���C���]���Ƃ���D
void B4P_Bearing::init_position(
	double balldia,			// �ʊ�{�a 
	int	   ballnum,			// �ʌ�
	double ballpcd,			// ��pcd
	double cos_alp0,		// �����ڐG�p�̗]��[-] 
	double omegair,			// ���։�]���x[rad/s]
	double omegaor,			// �O�։�]���x[rad/s]
	const double*rmg0		// �ێ���d�S����ɕێ���􉽒��S�֌������x�N�g��
) {
	// �����N�H�[�^�j�I���͑S�� [1,0,0,0] �Ƃ���D
	Quaterniond q0 = Quaterniond::Identity();

	// �ʒu�̎Z�o�D�ʂ�PCD��ɕ��Ԃ��̂Ƃ���D0�Ԗڂ̋ʂ�12���i+z���j�ɂ�����̂Ƃ���D
	double x; ArrayXd th(ballnum + 1), y(ballnum + 1), z(ballnum + 1);
	th = ArrayXd::LinSpaced(ballnum + 1, 0, 2 * Numeric::pi);
	x = 0;
	y = 0.5 * ballpcd *-th.sin();
	z = 0.5 * ballpcd * th.cos();

	// �������ʊp��ݒ肷��D
	this->azimuth0 = new double[ballnum];
	for (int i = 0; i < ballnum; i++)
		this->azimuth0[i] = th[i];

	// ���x�̎Z�o�D�ێ���̗��_���]���Ƃ���D�e�N�j�J�����|�[�gp.248�D�ʂ͐ڐG�_�ԋ����𒼌a�ɂ������Ɖ��肷��D
	double nc, Dw, vx; ArrayXd vy(ballnum + 1), vz(ballnum + 1);
	Dw = balldia * 0.5 * cos_alp0;	//�u�ڐG�_�ԋ����𒼌a�ɂ������v�̒��a

	nc = (1 - Dw / ballpcd) * 0.5 * omegair
		+ (1 + Dw / ballpcd) * 0.5 * omegaor; // ���_���]�� [rad/s]
	vx = 0;
	vy = 0.5 * ballpcd * nc *-th.cos();
	vz = 0.5 * ballpcd * nc *-th.sin();

	// �p���x�̎Z�o�D���_���]���Ƃ���D�e�N�j�J�����|�[�gp.248�D�ʂ͐ڐG�_�ԋ������������Ɖ��肷��D
	double na, wx, wy, wz;
	na = (ballpcd / Dw - Dw / ballpcd) * 0.5 * (omegaor - omegair); // ���_���]�� [rad/s]
	wx = na;
	wy = 0;
	wz = 0;

	Vector3d x0, v0, w0;
	// �e�]���̂ɒl�����D
	for (int i = 0; i < this->Z; i++) {
		x0(0) = x;		x0(1) = y[i];	x0(2) = z[i];
		v0(0) = vx;		v0(1) = vy[i];	v0(2) = vz[i];
		w0(0) = wx;		w0(1) = wy;		w0(2) = wz;

		this->BL[i].set_param(x0, v0, q0, w0);
	}

	// ���O�ցE�ێ���̏����ψʁE���x�E�p����0�D
	x0 *= 0;
	v0 *= 0;
	w0 *= 0;

	w0(0) = omegair;	// �p���x��x�����������͒l�D����0�D
	this->IR.set_param(x0, v0, q0, w0);

	w0(0) = omegaor;	// �p���x��x�����������͒l�D����0�D
	this->OR.set_param(x0, v0, q0, w0);

	x0 = - Vector3d(rmg0);	// �����ʒu�̏d�S�͕΂����ꏊ�ɑ��݂���D�􉽒��S�Əd�S�̍��͓��͒l�ł��邽�߂�������D
	w0(0) = nc;				// �p���x��x�����������_���]���D����0�D
	this->CG->set_param(x0, v0, q0, w0);

	return;
}


// �O���ϐ����玲��e���ނ̕ψʁE�^���ʂ�����������D
void B4P_Bearing::set_y(const double*y) {
	Vector3d    xi(y[0], y[1], y[2]);
	Vector3d    vi(y[3], y[4], y[5]);
	Quaterniond qi(y[6], y[7], y[8], y[9]);
	Vector3d    wi(y[10], y[11], y[12]);
	this->IR.set_y(xi, vi, qi, wi);
	//cout  << endl;
	//cout << "*qi=**************************" << endl;
	//std::cout << IR.q.x() << std::endl;
	//std::cout << IR.q.y() << std::endl;
	//std::cout << IR.q.z() << std::endl;
	//std::cout << IR.q.w() << std::endl;
	//cout << "wi=**************************" << endl;
	//
	//std::cout << IR.w << std::endl;
	//cout << "ax=**************************" << endl;
	//std::cout << IR.get_ax() << std::endl;

	Vector3d    xc(y[13], y[14], y[15]);
	Vector3d    vc(y[16], y[17], y[18]);
	Quaterniond qc(y[19], y[20], y[21], y[22]);
	Vector3d    wc(y[23], y[24], y[25]);
	this->CG->set_y(xc, vc, qc, wc);

	Vector3d xb, vb, wb;
	Quaterniond qb;
	for (int i = 0; i < this->Z; i++) {
		int j = 13 * i + 26;
		xb = Vector3d(y[j + 0], y[j + 1], y[j + 2]);
		vb = Vector3d(y[j + 3], y[j + 4], y[j + 5]);
		qb = Quaterniond(y[j + 6], y[j + 7], y[j + 8], y[j + 9]);
		wb = Vector3d(y[j + 10], y[j + 11], y[j + 12]);
		this->BL[i].set_y(xb, vb, qb, wb);
	}
	return;
}

// ����e���ނ̕ψʁE�^���ʂ��O���ϐ��ɏo�͂���D
void B4P_Bearing::get_y(double*y) {

	Vector3d x, p, w; Quaterniond q;

	for (int i = 0; i < this->Z; i++) {
		this->BL[i].get_y(x, p, q, w);
		y[13 * i + 26] = x.x(); y[13 * i + 27] = x.y(); y[13 * i + 28] = x.z();
		y[13 * i + 29] = p.x(); y[13 * i + 30] = p.y(); y[13 * i + 31] = p.z();
		y[13 * i + 32] = q.w(); y[13 * i + 33] = q.x(); y[13 * i + 34] = q.y(); y[13 * i + 35] = q.z();
		y[13 * i + 36] = w.x(); y[13 * i + 37] = w.y(); y[13 * i + 38] = w.z();
	}
	this->IR.get_y(x, p, q, w);
	y[0] = x.x(); y[1] = x.y(); y[2] = x.z();
	y[3] = p.x(); y[4] = p.y(); y[5] = p.z();
	y[6] = q.w(); y[7] = q.x(); y[8] = q.y(); y[9] = q.z();
	y[10] = w.x(); y[11] = w.y(); y[12] = w.z();

	this->CG->get_y(x, p, q, w);
	y[13] = x.x(); y[14] = x.y(); y[15] = x.z();
	y[16] = p.x(); y[17] = p.y(); y[18] = p.z();
	y[19] = q.w(); y[20] = q.x(); y[21] = q.y(); y[22] = q.z();
	y[23] = w.x(); y[24] = w.y(); y[25] = w.z();

	return;
}

// ���݂̏�Ԃł̊e���ނɂ�����׏d���v�Z���C�����l��Ԃ����\�b�h�D
void B4P_Bearing::get_dydt(double*dydt) {

	Vector3d sumFo, sumTo, sumFi, sumTi, sumFc, sumTc;
	sumFo = sumTo = sumFi = sumTi = sumFc = sumTc = Vector3d::Zero();

	for (int i = 0; i < this->Z; i++) {

		// �e�{�[���̉׏d�Ɩ��C���v�Z�D
		Vector3d Fbi, Tbi, Fib, Tib; // Fbi : �{�[��(b)������(i)����󂯂�́D
		this->BIP[i].calc_force(Fbi, Tbi, Fib, Tib);

		Vector3d Fbo, Tbo, Fob, Tob; // Fbo : �{�[��(b)���O��(o)����󂯂�́D
		this->BOP[i].calc_force(Fbo, Tbo, Fob, Tob);

		Vector3d Fbc, Tbc, Fcb, Tcb; // Fbc : �{�[��(b)���ێ���(c)����󂯂�́D
		this->BCP[i].calc_force(Fbc, Tbc, Fcb, Tcb);

		// �{�[���ɂ�����׏d�̑��a���Z�o�D
		Vector3d sumFb = Fbi + Fbo + Fbc;
		Vector3d sumTb = Tbi + Tbo + Tbc;

		// �d�͉����x�����Z�D
		sumFb += this->BL[i].get_mg();

		double dydt_[13];
		this->BL[i].get_dydt(sumFb, sumTb, dydt_);

		// ���[�v���Ƃɖ߂�l�ɑ���D
		for (int j = 0; j < 13; j++)
			dydt[13 * i + 26 + j] = dydt_[j];

		// �e���ނɂ�����׏d��ώZ�D
		sumFo += Fob;	sumTo += Tob;
		sumFi += Fib;	sumTi += Tib;
		sumFc += Fcb;	sumTc += Tcb;

		// �������ݗp�̕ϐ��ɕۑ��D
		this->BL[i].set_FT(sumFb, sumTb);
	}

	// �d�͉����x�����Z�D
	sumFo += this->OR.get_mg();
	sumFi += this->IR.get_mg();
	sumFc += this->CG->get_mg();

	// ���ւ������͉׏d�����Z�D�i������O�ցE�ێ���ɂ������\��H�j
	sumFi += this->F_load;
	sumTi += this->T_load;

	double dyidt_[13];
	this->IR.get_dydt_(sumFi, sumTi, dyidt_);

	for (int j = 0; j < 13; j++)
		dydt[0 + j] = dyidt_[j];

	double dycdt_[13];
	this->CG->get_dydt(sumFc, sumTc, dycdt_);

	for (int j = 0; j < 13; j++)
		dydt[13 + j] = dycdt_[j];

	// �������ݗp�̕ϐ��ɕۑ��D
	this->OR.set_FT(sumFo, sumTo);
	this->IR.set_FT(sumFi, sumTi);
	this->CG->set_FT(sumFc, sumTc);

	return;
}

// �O���ɏ������ނׂ��ϐ��� Eigen �̔z��ɂ��ĊO�ɓn�����\�b�h�D
void B4P_Bearing::save(B4P_Out&OUT) {

	for (int i = 0; i < this->Z; i++)
		this->BL[i].save(OUT.BL[i].x, OUT.BL[i].v, OUT.BL[i].q, OUT.BL[i].w, OUT.BL[i].ax, OUT.BL[i].F, OUT.BL[i].T);

	this->OR.save(OUT.OR.x, OUT.OR.v, OUT.OR.q, OUT.OR.w, OUT.OR.ax, OUT.OR.F, OUT.OR.T);
	this->IR.save(OUT.IR.x, OUT.IR.v, OUT.IR.q, OUT.IR.w, OUT.IR.ax, OUT.IR.F, OUT.IR.T);
	this->CG->save(OUT.CG.x, OUT.CG.v, OUT.CG.q, OUT.CG.w, OUT.CG.ax, OUT.CG.F, OUT.CG.T);
	
	for (int i = 0; i < this->Z; i++) {
		this->BOP[i].write(OUT.BOP[i]);
		this->BIP[i].write(OUT.BIP[i]);
		this->BCP[i].save(OUT.BCP[i]);
	}
	return;
}

// �É��(1) �����v�Z�F����e���ނ̕ψʁE�^���ʂ��O���ϐ��ɏo��
void B4P_Bearing::get_Xstf(double*Xstf) {
	this->IR.get_Xstf(Xstf);

	for (int i = 0; i < this->Z; i++) {
		Vector3d x, v, w;
		Quaterniond q;
		this->BL[i].get_param(x, v, q, w);
		Vector3d thXZ = this->OR.ine_to_XZcoord(x);
		Xstf[i * 2 + 5] = thXZ[1] / Rigid::l;
		Xstf[i * 2 + 6] = thXZ[2] / Rigid::l;
	}
	return;
}


// �É��(1) �����v�Z�Fcalculator�N���X����̓��� X_stf ������ɑ��
void B4P_Bearing::set_Xstf(double* Xstf) {

	this->IR.set_Xstf(Xstf);

	// �ʂ� ��-�� ���W�n���犵�����W�n�ɕϊ�
	for (int i = 0; i < this->Z; i++) {
		int i2 = i * 2;
		Vector3d eta(this->azimuth0[i], Xstf[i2 + 5] * Rigid::l, Xstf[i2 + 6] * Rigid::l);
		this->BOP[i].set_eta_stf(eta);
	}

	return;
}
// �É��(1) �����v�Z�F��������׏d���v�Z���Ccalculator�N���X�� F(X) �̌`�ŏo��
void B4P_Bearing::get_Fstf(double* Fstf) {

	Vector3d sumFi, sumTi;
	sumFi << 0, 0, 0;	sumTi << 0, 0, 0;

	// ���͉׏d�Əd�͉����x�����Z�D
	sumFi += this->F_load + this->IR.get_mg();
	sumTi += this->T_load;

	for (int i = 0; i < this->Z; i++) {

		// �e�{�[���̉׏d���v�Z�D
		Vector3d Fbi, Fbo;				// Fbi : �{�[��(b)������(i)����󂯂�́i�ʋO�����W�n�j
		Vector3d Fib, Tib, Fob, Tob;	// Fib : ����(i)���{�[��(b)����󂯂�́i�������W�n�j
		this->BIP[i].get_Fstf(Fbi, Fib, Tib);
		this->BOP[i].get_Fstf(Fbo, Fob, Tob);


		// �O�ցE���ւ��ꂼ��̍a���p�f�ʂŌ����Ƃ��̉׏d���v�Z�D
		Vector2d Fbo_ = this->OR.to_cross_sction(Fbo, this->azimuth0[i]);
		Vector2d Fbi_ = this->IR.to_cross_sction(Fbi, this->azimuth0[i]);
		Vector2d Fb = Fbo_ + Fbi_;

		// �ア�o�l�i�˗́j�Dx^-1 �̌`�ɂ��邱�Ƃɂ���āC�O���Ɏ������₷������D
		double x = Vector2d(this->BL[i].x[1], this->BL[i].x[2]).norm();
		double k = 1e-20;
		Fstf[2 * i + 5] = Fb[0] - k / x;
		Fstf[2 * i + 6] = Fb[1];

		// �ڐG�׏d�v�Z�ŏo�Ă������֑��̔���p��ώZ���Ă����D
		sumFi += Fib;
		sumTi += Tib;
	}



	// �S�Ẵ{�[���̉׏d�̑��a��߂�l�Ƃ��ė^����D
	Fstf[0] = sumFi[0];
	Fstf[1] = sumFi[1];
	Fstf[2] = sumFi[2];

	// �g���N����֍��W�n�ɕϊ����C��1,2������߂�l�ɂ���D�i��0�����͉񂷕����Ȃ̂ōl�����Ȃ��j
	Vector3d Ti = this->IR.to_myvector(sumTi);
	Fstf[3] = Ti[1];
	Fstf[4] = Ti[2];
	return;
}

B4P_Bearing::B4P_Bearing(void) {
	this->BL = NULL;
	this->BOP = NULL;
	this->BIP = NULL;
	this->BCP = NULL;
	this->azimuth0 = NULL;
	return;
}

B4P_Bearing::~B4P_Bearing(void) {
	if (this->BL != NULL)
		delete[] this->BL;
	if (this->BOP != NULL)
		delete[] this->BOP;
	if (this->BIP != NULL)
		delete[] this->BIP;
	if (this->BCP != NULL)
		delete[] this->BCP;
	if (this->azimuth0 != NULL)
		delete[] this->azimuth0;
	return;
}































//this->nX_stf =  2 * FI.ballnum + 5;
//this->nX_frc =  4 * FI.ballnum;
//this->nX_stf_frc =  6 * FI.ballnum + 5;


//// YZ�����O���׏d�̒P�ʃx�N�g�������߂�i1:Y�����C2:Z�����j
//Vector3d B4P_Bearing::get_Fyzload_dir() {
//	// YZ���ʏ�̗͂̑傫��
//	double F_str = sqrt(this->F_load[1]*this->F_load[1]+this->F_load[2]*this->F_load[2]);
//	if (F_str == 0) return Vector3d(0, 0, 0);
//	// �P�ʃx�N�g���̓��o
//	return Vector3d(0, this->F_load[1]/F_str, this->F_load[2]/F_str);
//}
//
//// X�����O���׏d�̌��������߂�i-1:�������C1:�������C0:X�׏d���j
//int B4P_Bearing::get_Fxload_dir() {
//
//	if (this->F_load[0] > 0) return 1;
//	else if (this->F_load[0] < 0) return -1;
//	else return 0;
//}
//
//// �ߎ���3.YZ�����O���׏d���[�����g�̕����x�N�g�������߂�i1:My�����C2:Mz�����j
//Vector3d B4P_Bearing::get_nyzload_dir() {
//	// YZ���ʏ�̗͂̑傫��
//	double N_str = sqrt(this->T_load[1]*this->T_load[1]+this->T_load[2]*this->T_load[2]);
//	if (N_str == 0) return Vector3d(0, 0, 0);
//	// �P�ʃx�N�g���̓��o
//	return Vector3d(0, this->T_load[1]/N_str, this->T_load[2]/N_str);
//}
//
//
//
//// �ߎ���3. �������̎p�������
//void B4P_Bearing::set_q_ir(Vector3d ax) {
//	this->IR.set_ax(ax);
//	return;
//}

//
//
//// ���݂̏�Ԃł̓��ցE�ʂP�ɂ�����׏d���Z�o���郁�\�b�h�D�i�S�Ċ������W�n�DF[0~2]�F���։׏d�CF[3~4]�F���փg���N�CF[5~6]�F�ʉ׏d�D�j
//void B4P_Bearing::get_F_stf_ball(double*F, int i) {
//
//	// �e�{�[���̉׏d���v�Z�D
//	Vector3d Fbi, Fib, Tib; // Fbi : �{�[��(b)������(i)����󂯂�́D
//	this->BIP[i].calc_force_stf(Fbi, Fib, Tib);
//	Vector3d Fbo, Fob, Tob; // Fbo : �{�[��(b)���O��(o)����󂯂�́D
//	this->BOP[i].calc_force_stf(Fbo, Fob, Tob);
//
//	// �d�͉����x�����Z�D���������O�֊�Ƃ������ƂŁC�O�։׏d�ɉ��Z���Ă���D
//	// ����ɉ��S�͂����Z�D
//	Fbo += this->BL[i].get_mg() + this->BOP[i].get_centf();
//
//	// �O�ցE���ւ��ꂼ��̍a���p�f�ʂŌ����Ƃ��̉׏d���v�Z�D
//	Vector2d Fbo_ = this->OR.to_cross_sction(Fbo, this->azimuth[i]);
//	Vector2d Fbi_ = this->IR.to_cross_sction(Fbi, this->azimuth[i]);
//	Vector2d Fb   = Fbo_ + Fbi_;
//
//
//
//
//	// �{�[���ɂ�����׏d�̑��a���Z�o�D
//	F[5] = Fb[0];
//	F[6] = Fb[1];
//
//	// �S�Ẵ{�[���̉׏d�̑��a��߂�l�Ƃ��ė^����D
//	F[0] = Fib[0];
//	F[1] = Fib[1];
//	F[2] = Fib[2];
//
//	// �g���N����֍��W�n�ɕϊ����C��1,2������߂�l�ɂ���D�i��0�����͉񂷕����Ȃ̂ōl�����Ȃ��j
//	Vector3d Ti = this->IR.to_myvector(Tib);
//	F[3] = Ti[1];
//	F[4] = Ti[2];
//}
//
//
//// �O���ϐ����玲��e���ނ̕ψʁE�^���ʂ�����������D�i���������v�Z�FSTEP1�p�D�j
//void B4P_Bearing::set_param_stf(const double*y_stf) {
//	//���ւ̊e�ϐ�������������
//	Vector3d    xi(y_stf[0], y_stf[1], y_stf[2]);
//	Vector3d    v0(0.0, 0.0, 0.0);
//	Quaterniond qi(1.0, 0.0, y_stf[3], y_stf[4]);
//	//Vector3d    w0( 1000.0,  0.0,  0.0);
//	Vector3d    w0(0.0, 0.0, 0.0);
//
//	this->IR.set_param(xi, v0, qi, w0);
//
//	for (int i = 0; i < this->Z; i++) {
//		Vector3d xb = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5], y_stf[2*i+6]));
//		this->BL[i].set_param_stf(xb);
//	}
//}
//
//
//// �ψʂ� y_stf �ł���Ƃ��� Jacobian ���Z�o���郁�\�b�h�D
//void B4P_Bearing::get_Jacobian_stf(double*y_stf, double dx, double dq, double*Jacobian) {
//
//	int N = this->nX_stf;
//
//	for (int i = 0; i < N * N; i++)
//		Jacobian[i] = 0.0;
//
//	double F_Xp[7], F_Xm[7], F_Zp[7], F_Zm[7];
//	for (int i = 0; i < this->Z; i++) {
//		Vector3d x = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5], y_stf[2*i+6]));
//		this->BL[i].set_param_stf(x);
//	}
//
//
//	for (int i = 0; i < this->Z; i++) {
//		// X�����ɔ����ψ�+dx��ǉ�
//		Vector3d Xp = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5]+dx, y_stf[2*i+6]));
//		this->BL[i].set_param_stf(Xp);
//		this->get_F_stf_ball(F_Xp, i);
//
//		// X�����ɔ����ψ�-dx��ǉ�
//		Vector3d Xm = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5]-dx, y_stf[2*i+6]));
//		this->BL[i].set_param_stf(Xm);
//		this->get_F_stf_ball(F_Xm, i);
//
//		// X�ɂ��ĕΔ������Ƃ�D(F[0~2]�F���։׏d�CF[3~4]�F���փg���N�CF[5~6]�F�ʉ׏d�C�Ƃ������Ƃ��l���ɓ���Čv�Z����j
//		//for (int j = 0; j < 5; j++){
//		int j;
//		j = 0;	Jacobian[N * j + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		j = 1;	Jacobian[N * j + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		j = 2;	Jacobian[N * j + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		j = 3;	Jacobian[N * j + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		j = 4;	Jacobian[N * j + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		j = 5;	Jacobian[N * (2 * i + 5) + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		j = 6;	Jacobian[N * (2 * i + 6) + 2 * i + 5] = (F_Xp[j] - F_Xm[j]) / (2 * dx);
//		//}
//
//		// Z�����ɔ����ψ�+dx��ǉ�
//		Vector3d Zp = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5], y_stf[2*i+6]+dx));
//		this->BL[i].set_param_stf(Zp);
//		this->get_F_stf_ball(F_Zp, i);
//
//		// Z�����ɔ����ψ�-dx��ǉ�
//		Vector3d Zm = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5], y_stf[2*i+6]-dx));
//		this->BL[i].set_param_stf(Zm);
//		this->get_F_stf_ball(F_Zm, i);
//
//		// Z�ɂ��ĕΔ������Ƃ�D
//		//int j;
//		j = 0;	Jacobian[N * j + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//		j = 1;	Jacobian[N * j + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//		j = 2;	Jacobian[N * j + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//		j = 3;	Jacobian[N * j + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//		j = 4;	Jacobian[N * j + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//		j = 5;	Jacobian[N * (2 * i + 5) + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//		j = 6;	Jacobian[N * (2 * i + 6) + 2 * i + 6] = (F_Zp[j] - F_Zm[j]) / (2 * dx);
//
//
//		// �{�[�������̈ʒu�ɖ߂��D
//		Vector3d x = this->OR.XZ_to_inecoord(Vector3d(this->azimuth[i], y_stf[2*i+5], y_stf[2*i+6]));
//		this->BL[i].set_param_stf(x);
//	}
//
//	Vector3d    xi(y_stf[0], y_stf[1], y_stf[2]);
//	Vector3d    v0(0.0, 0.0, 0.0);
//	Quaterniond qi(1.0, 0.0, y_stf[3], y_stf[4]);
//	Vector3d    w0(0.0, 0.0, 0.0);
//	double*F_p = new double[N];
//	double*F_m = new double[N];
//
//	for (int i = 0; i < 3; i++) {
//		xi[i] = y_stf[i] + dx;
//		this->IR.set_param(xi, v0, qi, w0);
//		this->get_F_stf(F_p);
//
//		xi[i] = y_stf[i] - dx;
//		this->IR.set_param(xi, v0, qi, w0);
//		this->get_F_stf(F_m);
//
//		xi[i] = y_stf[i];
//		this->IR.set_param(xi, v0, qi, w0);
//
//		for (int j = 0; j < N; j++) {
//			Jacobian[N * i + j] = (F_p[j] - F_m[j]) / (2 * dx);
//		}
//	}
//
//	// �N�H�[�^�j�I����y,z���ꂼ��v�Z�D
//	qi.y() = y_stf[3] + dx;
//	this->IR.set_param(xi, v0, qi, w0);
//	this->get_F_stf(F_p);
//
//	qi.y() = y_stf[3] - dx;
//	this->IR.set_param(xi, v0, qi, w0);
//	this->get_F_stf(F_m);
//
//	for (int j = 0; j < N; j++) {
//		Jacobian[N * 3 + j] = (F_p[j] - F_m[j]) / (2 * dx);
//		//cout << N * 3 + j << ":" <<
//		//	Jacobian[N * 3 + j] << endl;
//	}
//	qi.y() = y_stf[3];
//	qi.z() = y_stf[4] + dx;
//	this->IR.set_param(xi, v0, qi, w0);
//	this->get_F_stf(F_p);
//
//	qi.z() = y_stf[4] - dx;
//	this->IR.set_param(xi, v0, qi, w0);
//	this->get_F_stf(F_m);
//
//	for (int j = 0; j < N; j++) {
//		Jacobian[N * 4 + j] = (F_p[j] - F_m[j]) / (2 * dx);
//		//cout << N * 4 + j << ":" <<
//		//	Jacobian[N * 4 + j] << endl;
//	}
//	// ���ɖ߂�
//	this->set_param_stf(y_stf);
//
//}
//
//



//// �O���ϐ����玲��e���ނ̕ψʁE�^���ʂ�����������D�i���C�����v�Z�FSTEP2�p�D�j
//void B4P_Bearing::set_param_frc(const double*y_frc) {
//
//	for (int i = 0; i < this->Z; i++) {
//		int j = 4 * i;
//		Vector3d vb = this->OR.Y_to_inevelocity(y_frc[j+0], this->azimuth[i]);
//		Vector3d wb = Vector3d(y_frc[j+1], y_frc[j+2], y_frc[j+3]);
//		this->BL[i].set_param_frc(vb, wb);
//	}
//}
//
//// �ێ�����]�����擾���郁�\�b�h�D
//double B4P_Bearing::get_wCage(void) {
//	return -1e5;
//}
//
//// ����e���ނ̕ψʁE�^���ʂ��O���ϐ��ɏo�͂���D�i���C�����v�Z�FSTEP2�p�D�j
//void B4P_Bearing::get_param_frc(double*y_frc) {
//	Vector3d vb, wb;
//	for (int i = 0; i < this->Z; i++) {
//		int j = 4 * i;
//		this->BL[i].get_param_frc(vb, wb);
//		double v = this->OR.ine_to_Yvelocity(vb, this->azimuth[i]);
//		y_frc[j+0] = v; y_frc[j+1] = wb[0]; y_frc[j+2] = wb[1]; y_frc[j+3] = wb[2];
//	}
//}
	/*��input�t�@�C������̓ǂݍ���*/
	// �{�[�����������o���C�z��̓��I�m�ہD
	/*
	int Z     = FI.Z;
	this->Z   = Z;
	this->BL  = new Ball[Z];
	this->BOP = new B4P_BallRingPair[Z];
	this->BIP = new B4P_BallRingPair[Z];
	this->BCP = new B4P_BallCagePair[Z];
	*/
	/*�Vinput�t�@�C������̓ǂݍ���*/
