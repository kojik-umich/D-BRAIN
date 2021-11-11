#include "bal_Cage.h"

// ���ەێ���N���X�D
void bal_Cage::Nothing(void) {
}

// ���^�ێ���N���X�D
void bal_SnapCage::init(const B4P_In&FI) {

	this->ro = FI.Snap.dout * 0.5;
	this->ri = FI.Snap.din * 0.5;
	this->m = FI.Cage.m;
	this->pcd = this->ro + this->ri + FI.Snap.jc;
	this->m_inv = 1.0 / this->m;
	for (int i = 0; i < 3; i++)
		this->xc[i] = FI.Cage.rmg0[i];
	this->Z = FI.ballnum;
	this->I = Vector3d(FI.Cage.Ix, FI.Cage.Iyz, FI.Cage.Iyz);
	this->I_inv = Vector3d(1 / FI.Cage.Ix, 1 / FI.Cage.Iyz, 1 / FI.Cage.Iyz);

	double x; ArrayXd th(this->Z + 1), y(this->Z + 1), z(this->Z + 1);
	th = ArrayXd::LinSpaced(this->Z + 1, 0, 2 * Numeric::pi);
	x = 0;
	y = 0.5 * this->pcd *-th.sin();
	z = 0.5 * this->pcd * th.cos();

	// �|�P�b�g�z��𓮓I�m�ہD
	this->PK = new Pocket[FI.ballnum];

	// �S�|�P�b�g���������D
	for (int i = 0; i < Z; i++) {
		Vector3d ax = this->ax0;		// ������
		this->PK[i].R = FI.Snap.R;
		this->PK[i].x = Vector3d(x, y[i], z[i]);
		this->PK[i].er = this->PK[i].x.normalized();
		this->PK[i].eth = this->PK[i].er.cross(ax);
		this->PK[i].x_norm = this->PK[i].x.norm();

		this->PK[i].ropen = FI.Snap.ropen;
		this->PK[i].kface = FI.Snap.Kface;			// �ʐڐG����[N/m]
		this->PK[i].kci = FI.Snap.Kcornerin;		// �p�ڐG����[N/m]
		this->PK[i].kco = FI.Snap.Kcornerout;		// �p�ڐG����[N/m]
		this->PK[i].kopen = FI.Snap.Kopen;			// �J��������[N/m]
		this->PK[i].kedgei = FI.Snap.Kedgein;			// �G�b�W����[N/m]	
		this->PK[i].kedgeo = FI.Snap.Kedgeout;			// �G�b�W����[N/m]

		// �J�����J���p��
		this->PK[i].cos_alp = sqrt(1.0 - FI.Snap.ropen*FI.Snap.ropen / this->PK[i].R / this->PK[i].R);
		this->PK[i].h0 = this->PK[i].R * this->PK[i].cos_alp;

		// �J�����J���p��
		double cosi = Numeric::Law_of_cosines(this->PK[i].ropen, ri, this->PK[i].x_norm);
		double coso = Numeric::Law_of_cosines(this->PK[i].ropen, ro, this->PK[i].x_norm);
		double sini = sqrt(1.0 - cosi * cosi);
		double sino = sqrt(1.0 - coso * coso);
		this->PK[i].cos_gammai = (ri * cosi - this->PK[i].x_norm) / sqrt(pow(ri * cosi - this->PK[i].x_norm, 2) + pow(ri * sini, 2));
		this->PK[i].cos_gammao = (ro * coso - this->PK[i].x_norm) / sqrt(pow(ro * coso - this->PK[i].x_norm, 2) + pow(ro * sino, 2));

		// ���������̌v�Z
		double cosim = Numeric::Law_of_cosines(this->PK[i].R, ri, this->PK[i].x_norm);
		this->PK[i].zimax = this->ri;
		this->PK[i].zimin = this->ri * cosim;

		double cosom = Numeric::Law_of_cosines(this->PK[i].R, ro, this->PK[i].x_norm);
		this->PK[i].zomax = this->ro;
		this->PK[i].zomin = this->ro * cosom;

		this->PK[i].ziave = (this->PK[i].zimax + this->PK[i].zimin) * 0.5;
		this->PK[i].ziamp = (this->PK[i].zimax - this->PK[i].zimin) * 0.5;

		this->PK[i].zoave = (this->PK[i].zomax + this->PK[i].zomin) * 0.5;
		this->PK[i].zoamp = (this->PK[i].zomax - this->PK[i].zomin) * 0.5;

		// �R�[�i�[�̍��W�̐ݒ�D
		this->PK[i].corni0 = this->PK[i].h0 * ax - ri * sini * this->PK[i].eth + ri * cosi * this->PK[i].er;
		this->PK[i].corni1 = this->PK[i].h0 * ax + ri * sini * this->PK[i].eth + ri * cosi * this->PK[i].er;
		this->PK[i].corno0 = this->PK[i].h0 * ax - ro * sino * this->PK[i].eth + ro * coso * this->PK[i].er;
		this->PK[i].corno1 = this->PK[i].h0 * ax + ro * sino * this->PK[i].eth + ro * coso * this->PK[i].er;


	}
	return;
}

// �|�P�b�g���S����Ƃ����ʕ����x�N�g������ڐG�p�^�[���𔻒�i�ڐG�̗L���͔��肵�Ȃ��j
bal_SnapCage::ContactPattern bal_SnapCage::how_Contact
(							// out:	[-]:	�ڐG�p�^�[���D
	const Vector3d&BL_x_,	// in:	[m]:	�{�[���̈ʒu�i�ێ���􉽒��S���W�n�j�D
	int np,					// in:	[-]:	���̃{�[���̑Ή�����|�P�b�g�ԍ��D
	Vector3d&edir,			// out: [-]:	�|�P�b�g���S����ʒ��S�֌������P�ʃx�N�g���D
	double&zo_th,			// out: [m]:	�ێ��풆�S����|�P�b�g���S�����Ō����C���݂̕��ʊp�ɂ�����ێ���O�G�b�W�̍����D
	double&zi_th, 			// out: [m]:	�ێ��풆�S����|�P�b�g���S�����Ō����C���݂̕��ʊp�ɂ�����ێ�����G�b�W�̍����D
	double&cos_th, 			// out: [-]:	��������12���Ɍ����Ă��ۂ̕��ʊp��cosine�D
	double&sin_th, 			// out: [-]:	��������12���Ɍ����Ă��ۂ̕��ʊp��sine�D
	double&cos_gamma, 		// out: [-]:	����������݂��|�P�b�g�𒆐S�Ƃ����ʂ̈ʑ��p��cosine�D
	double&sin_gamma		// out: [-]:	����������݂��|�P�b�g�𒆐S�Ƃ����ʂ̈ʑ��p��sine�D
) {
	// �O�����Ƃ��ă|�P�b�g���S����Ƃ����ʕ����x�N�g�� edir ���v�Z�D�i�ێ���􉽒��S���W�n�j
	Vector3d pb = BL_x_ - this->PK[np].x;
	double pb_norm = pb.norm();

	// 0���Z�̏ꍇ�͗�O�����i�|�P�b�g���S�ɋʂ�����ꍇ�j
	if (pb_norm == 0.0)
		return exception;
	edir = pb / pb_norm;

	// (1) �ʕ��ʊp���𓱏o�i����������Ƃ����Ƃ��̕��ʊp�C�J���������蔻��Ɏg�p�j
	// ���ʊp���̗]�����v�Z�D
	double   cos_alp = this->ax0.dot(edir);	// ���F��=[0,180]�Ȃ̂�cos�ň�ӂɌ��܂�D

	// (2) �ʕ��ʊp���𓱏o�i"����������݂�"�a��������Ƃ����Ƃ��̕��ʊp�C�p�����蔻��Ɏg�p�j
	// �O�ς�p���Ča������ ��+90�� ���Ȃ��x�N�g�����v�Z
	Vector3d edirax = edir.cross(this->ax0);
	double edirax_norm = edirax.norm();

	// 0���Z�̏ꍇ�͗�O�����i���ƃ|�P�b�g���S-�{�[�����S�����s�j
	if (edirax_norm == 0.0)
		return ContactPattern::exception;

	Vector3d azi = edirax / edirax_norm;
	cos_gamma = this->PK[np].eth.dot(azi);
	sin_gamma = -this->PK[np].er.dot(azi);

	// (3) ���������̌v�Z�i�G�b�W������Ɏg�p�j
	// �ێ���􉽒��S�E�|�P�b�g���S�E�{�[�����S�̒���O�p�`�̖@���x�N�g���D�i�ێ���􉽒��S���W�n�j
	Vector3d eredir = this->PK[np].er.cross(edir);
	double   eredir_norm = eredir.norm();

	// 0���Z�̏ꍇ�͗�O�����i�ێ���􉽒��S�E�|�P�b�g���S�E�{�[�����S�����j
	if (eredir_norm == 0.0)
		return exception;

	Vector3d en = eredir / eredir_norm;

	// ���ʊp �� �� sine.�i12�����ێ���􉽒��S���W�n��[1,0,0]�����C3�����������D�j
	sin_th = -this->ax0.dot(en);
	double cos_2th = 1.0 - 2 * sin_th * sin_th;
	cos_th = this->PK[np].eth.dot(en);

	// �ߎ����ɂ���ĕ��ʂɂ�鐂������������D�i�ʌa�ɑ΂��Č덷 0.01% ���x���m�F�ς݁D�j
	zi_th = this->PK[np].ziave + this->PK[np].ziamp * cos_2th;
	zo_th = this->PK[np].zoave + this->PK[np].zoamp * cos_2th;

	// �ڐG������̐������� L �̎Z�o�D
	double cos_beta = this->PK[np].er.dot(edir);
	double L = this->PK[np].x_norm + this->PK[np].R * cos_beta;

	// �Ō�ɏꍇ���������{�D
	if (cos_alp < this->PK[np].cos_alp) {
		if (L < zi_th)
			return edgei;
		else if (L > zo_th)
			return edgeo;
		else
			return face;
	}
	if (cos_gamma < this->PK[np].cos_gammai)
		return corneri;
	if (cos_gamma > this->PK[np].cos_gammao)
		return cornero;

	return aperture;
}

// �ڐG�p�^�[���ɉ����ĐڐG�_�����v�Z���郁�\�b�h�D�ێ���͍��̂Ɖ���D�i�ێ���􉽒��S���W�n�j
int bal_SnapCage::where_Contact
(							// out:	[-]:	�ڐG�_���D
	ContactPattern c,		// in:	[-]:	�ڐG�p�^�[���D
	int np,					// in:	[-]:	���̃{�[���̑Ή�����|�P�b�g�ԍ��D
	const Vector3d&dir, 	// in:	[-]:	�|�P�b�g���S����ʒ��S�֌������P�ʃx�N�g���D
	double zo_th, 			// in:	[m]:	�ێ��풆�S����|�P�b�g���S�����Ō����C���݂̕��ʊp�ɂ�����ێ���O�G�b�W�̍����D
	double zi_th, 			// in:	[m]:	�ێ��풆�S����|�P�b�g���S�����Ō����C���݂̕��ʊp�ɂ�����ێ�����G�b�W�̍����D
	double cos_th, 			// in:	[-]:	��������12���Ɍ����Ă��ۂ̕��ʊp��cosine�D
	double sin_th, 			// in:	[-]:	��������12���Ɍ����Ă��ۂ̕��ʊp��sine�D
	double cos_gamma, 		// in:	[-]:	����������݂��|�P�b�g�𒆐S�Ƃ����ʂ̈ʑ��p��cosine�D
	double sin_gamma, 		// in:	[-]:	����������݂��|�P�b�g�𒆐S�Ƃ����ʂ̈ʑ��p��sine�D
	Vector3d * x_,			// out:	[m]:	�ڐG�_�ʒu�i�ێ���􉽒��S���W�n�j�D�ڐG�_���������|�C���^�ŕ�
	double * k	   			// out:	[N/m]:	�e�ڐG�_�ɂ����鍄���D�ڐG�_���������|�C���^�ŕԂ��D
) {
	switch (c) {
	case face: {
		x_[0] = this->PK[np].x + dir * this->PK[np].R;
		k[0] = this->PK[np].kface;
		return 1;
	}
	case edgeo: {
		double R_cs = sqrt(Numeric::Square(this->PK[np].R) - Numeric::Square(zo_th - this->PK[np].x_norm));
		x_[0] = this->ax0    * cos_th * R_cs
			+ this->PK[np].eth * sin_th * R_cs
			+ this->PK[np].er * zo_th;
		k[0] = this->PK[np].kedgeo;
		x_[1] = - this->ax0    * cos_th * R_cs
			- this->PK[np].eth * sin_th * R_cs
			+ this->PK[np].er * zo_th;
		k[1] = this->PK[np].kedgeo;
		return 2;
	}
	case edgei: {
		double R_cs = sqrt(Numeric::Square(this->PK[np].R) - Numeric::Square(zi_th - this->PK[np].x_norm));
		x_[0] = this->ax0    * cos_th * R_cs
			+ this->PK[np].eth * sin_th * R_cs
			+ this->PK[np].er * zi_th;
		k[0] = this->PK[np].kedgei;
		x_[1] = - this->ax0    * cos_th * R_cs
			- this->PK[np].eth * sin_th * R_cs
			+ this->PK[np].er * zi_th;
		k[1] = this->PK[np].kedgei;
		return 2;
	}
	case aperture: {
		x_[0] = this->PK[np].x
			+ this->ax0      * this->PK[np].h0
			+ this->PK[np].eth * sin_gamma * this->PK[np].ropen
			+ this->PK[np].er * cos_gamma * this->PK[np].ropen;
		k[0] = this->PK[np].kopen;
		x_[1] = this->PK[np].x
			+ this->ax0      * this->PK[np].h0
			- this->PK[np].eth * sin_gamma * this->PK[np].ropen
			- this->PK[np].er * cos_gamma * this->PK[np].ropen;
		k[1] = this->PK[np].kopen;
		return 2;
	}
	case cornero: {
		x_[0] = this->PK[np].corno0;
		x_[1] = this->PK[np].corno1;
		k[0] = this->PK[np].kco;
		k[1] = this->PK[np].kco;
		return 2;
	}
	case corneri: {
		x_[0] = this->PK[np].corni0;
		x_[1] = this->PK[np].corni1;
		k[0] = this->PK[np].kci;
		k[1] = this->PK[np].kci;
		return 2;
	}
	case exception:
		return 0;

	default:
		return 0;
	}
}


// �{�[���̈ʒu�Ɣ��a�ɂ���ĐڐG�_�ʒu�ƍ�����Ԃ����\�b�h�D���̏������ƂɊO���ŉ׏d���C�v�Z���s���D
int bal_SnapCage::get_ContactPoint(
	const Vector3d&BL_x,	// in:	[m]:	�{�[���̈ʒu�i�������W�n�j�D
	double         BL_r,	// in:	[m]:	�{�[���̔��a�D 
	int            np, 		// in:	[-]:	���̃{�[���̑Ή�����|�P�b�g�ԍ��D
	Vector3d      *x, 		// out:	[m]:	�ڐG�_�ʒu�i�������W�n�j�D�ڐG�_���������|�C���^�ŕԂ��D
	double        *k, 		// out:	[N/m]:	�e�ڐG�_�ɂ����鍄���D�ڐG�_���������|�C���^�ŕԂ��D
	int			  *ptt		// out: [-]:	�ڐG�p�^�[��(int�^�ɕϊ����ďo��)
) {
	// �{�[���ʒu�����߂�D�i�ێ���􉽒��S���W�n�j
	Vector3d BL_x_ = this->to_mycoord(BL_x) - this->xc;

	Vector3d dir = Vector3d::Zero();
	double zo_th, zi_th, cos_th, sin_th, cos_gamma, sin_gamma;
	zo_th = zi_th = cos_th = sin_th = cos_gamma = sin_gamma = 0.0;

	ContactPattern c = this->how_Contact(BL_x_, np, dir, zo_th, zi_th, cos_th, sin_th, cos_gamma, sin_gamma);

	Vector3d x_[_MAX_CONTACT_];	// �ڐG�_�ʒu�̔z��D�i�ێ���􉽒��S���W�n�j
	int nc = this->where_Contact(c, np, dir, zo_th, zi_th, cos_th, sin_th, cos_gamma, sin_gamma, x_, k);

	// �ێ���􉽒��S���W�n �� �ێ�����W�n �� �������W�n�ɕϊ��D
	for (int i = 0; i < nc; i++) {
		x[i] = this->to_inecoord(x_[i] + this->xc);
		ptt[i] = static_cast<int>(c);
	}
	return nc;
}

bal_SnapCage::bal_SnapCage() {
	this->PK = NULL;
	return;
}

bal_SnapCage::~bal_SnapCage() {
	if (this->PK != NULL)
		delete[] this->PK;
	return;
}



//���^�ێ���d�ʎ����v�Z
//double bal_SnapCage::calc_snapmass(double Snap_den_, double Snap_rout_, double Snap_rin_, double Snap_h_, double ballpcd_, double balldia_, double clr_, double ballnum_){
//	double rp3 = 0.5*ballpcd_+0.5*balldia_+0.25*clr_;//�O�֔��a
//	double mc;
//	/* snap cage 1.4465 is modification factor */
//	mc=1.4465*Snap_den_*Numeric::pi*((Snap_rout_*Snap_rout_-Snap_rin_*Snap_rin_)*Snap_h_
//		-(double)ballnum_*(Snap_rout_-Snap_rin_)*rp3*rp3
//		*(0.5+0.5*sin(2.0*asin((Snap_h_)/rp3))+asin((Snap_h_)/rp3)));
//	return mc;
//}

