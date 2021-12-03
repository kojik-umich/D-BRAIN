/*******************************************************************************
!								"BS_BallScrew.cpp"
!													2020/04/15	[Core-T]	���
!
! �{�[���˂��I�u�W�F�N�g�D���Ō����Ƃ���̎���I�u�W�F�N�g�ɑ����D
!
!*******************************************************************************/

#include "BS_BallScrew.h"

// �{�[���˂��\���v�f�𓮓I�m�ۂ���D��ԍŏ��ɋN�����Ăق������\�b�h�D
void BS_BallScrew::allocate(const BS_In&IN) {

	if (IN.nut.size() == 1)
		this->NT = new BS_SingleNut();
	else
		this->NT = new BS_DoubleNut();

	this->nCC = IN.circuit.size();
	this->CC = new BS_Circuit[this->nCC];
	for (size_t i = 0; i < IN.circuit.size(); i++)
		this->CC[i].allocate(IN.circuit[i].ball.size());

	this->nP = 0;
	for (size_t i = 0; i < IN.circuit.size(); i++)
		this->nP += IN.circuit[i].ball.size();

	this->BNP = new BS_BallNutPair[this->nP];
	this->BSP = new BS_BallShaftPair[this->nP];

	this->NT->allocate(IN.nut);
	this->ST.allocate(IN.shaft);

	return;
}

// �{�[���ƑΉ����镔�ނ̔ԍ���R�t����D
void BS_BallScrew::link(const BS_In&IN) {

	int iter = 0;

	for (int i = 0; i < this->nCC; i++) {

		this->CC[i].link(
			&this->NT->CY[IN.circuit[i].inut],
			IN.circuit[i].is);

		for (int j = 0; j < this->CC[i].nBL; j++) {

			this->BNP[iter].link(
				&this->CC[i].BL[j],
				&this->NT->CY[IN.circuit[i].inut],
				IN.circuit[i].is);


			this->BSP[iter].link(
				&this->CC[i].BL[j],
				&this->ST.CY[0],
				IN.circuit[i].is);
			iter++;
		}
	}
	return;
}

// ��H�ɂ���e�ʂ̍��W�����������郁�\�b�h�D�Ƃ肠����0�n�܂�Ŏ����Ă݂�D
void BS_BallScrew::init_position(
	double wn,	// in:	[rad/s]		���_���]���D
	double ws	// in:	[rad/s]		���_���]���D
) {
	for (int i = 0; i < this->nCC; i++) {

		Vector3d dx = this->CC[i].CY->x;

		Vector3d*x0s = new Vector3d[this->CC[i].nBL];
		this->CC[i].get_x0s(x0s);

		for (int j = 0; j < this->CC[i].nBL; j++)
			this->CC[i].BL[j].x = x0s[j] + 0.5 * dx;

		double nd = this->CC[i].get_nd();
		double pcd = this->CC[i].get_r() * 2;

		for (int j = 0; j < this->CC[i].nBL; j++) {
			Vector3d t = this->CC->get_t(j);
			double D = this->CC[i].BL[j].r * sqrt(2.0);	// �ڐG�p45�x������D
			double w = 0.5 * (1.0 - D / pcd) * ws + 0.5 * (1.0 + D / pcd) * wn;	// ���]�p���x
			Vector3d v = w * nd * t;
			Vector3d w_ = (pcd / D - D / pcd) * 0.5 * (wn - ws) * Vector3d(1.0, 0.0, 0.0);	// ���]�p���x

			this->CC[i].BL[j].v = v;
			this->CC[i].BL[j].w = w_;
			this->CC[i].BL[j].q.setIdentity();
		}
	}
	return;
}
// �׏d�̏������D�܂��C�׏d���烂�[�����g���v�Z���C���̓��[�����g�ɉ��Z�D
void BS_BallScrew::init_Load(const BS_In&IN) {

	Ball Dummy;
	Dummy.x = Vector3d::Zero();
	this->LD.F = Vector3d::Zero();
	this->LD.T = Vector3d::Zero();

	for (size_t i = 0; i < IN.load.size(); i++) {

		this->LD.F += Vector3d(IN.load[i].F);
		this->LD.T += Dummy.calc_Torque(
			Vector3d(IN.load[i].x),
			Vector3d(IN.load[i].F))
			+ Vector3d(IN.load[i].T);
	}
	return;
}

// ���͂���{�[���˂��̏��������s���D
void BS_BallScrew::init(const BS_In&IN, double v0, double w0, double wn) {

	this->link(IN);

	for (int i = 0; i < this->nCC; i++)
		this->CC[i].init(IN.circuit[i]);

	this->NT->init(IN.nut, wn);
	this->ST.init(IN.shaft, IN.bound.v_const, IN.bound.w_const, IN.bound.tan_thy, IN.bound.tan_thz, v0, w0);

	for (int i = 0; i < this->nP; i++) {
		this->BNP[i].init(IN.BallNutPair[i], IN.tribology, IN.oil);
		this->BSP[i].init(IN.BallShaftPair[i], IN.tribology, IN.oil);
	}
	this->init_position(wn, w0);
	this->ST.init_pos(v0, w0);

	for (int i = 0; i < this->nP; i++) {
		Vector3d etan = this->BNP[i].get_eta();
		this->BNP[i].th0 = etan[0];
		Vector3d etas = this->BSP[i].get_eta();
		this->BSP[i].thp = etas[0];
	}
	this->init_Load(IN);

	for (int j = 0; j < this->nCC; j++) {
		this->CC[j].iSP = IN.circuit[j].is;
		this->CC[j].nBL = IN.circuit[j].ball.size();

		this->CC[j].th0 = IN.circuit[j].th0;
		this->CC[j].th1 = IN.circuit[j].th1;
	}
	return;
}

// �y�����l�v�Z�z�V���t�g�C�{�[���̏����ʒu���������茈�߂郁�\�b�h�D
void BS_BallScrew::preset_y0(double dx0, double dth0, double dx1) {

	this->preset_y0_F(dx0);
	this->preset_y0_T(dth0);
	this->preset_y0_x(dx1);

	return;
}

// ���͂ɂ��鏉���ʒu�ɃV���t�g��ݒ肷�郁�\�b�h�D
void BS_BallScrew::lock_y0(const double*x0, const double*ax0, double v0, double w0) {

	Vector3d x = Vector3d(x0);
	Vector3d ax = Vector3d(ax0);

	this->ST.x = x;
	this->ST.set_ax(ax);
	this->ST.v = v0 * ax;
	this->ST.w = w0 * ax;

	this->ST.set_dx();

	for (int i = 0; i < this->nCC; i++)
		for (int j = 0; j < this->CC[i].nBL; j++)
			this->CC[i].BL[j].x += Vector3d(x0) / 2;

	return;
}

// �y�����l�v�Z�z�{�[���S�̂ƃV���t�g���׏d�����ɓ������āC�������菉���ʒu�����߂郁�\�b�h�D
void BS_BallScrew::preset_y0_F(
	double dx0		// in:	[m]	���������ω��ʁD�ʏ�1nm�Ȃǌn�ɑ΂��ď\���������l��p���Ă��������D
) {
	// �׏d���\���������ꍇ�C���̃��\�b�h�͎g���Ȃ��D
	double F_norm = this->LD.F.norm();
	if (F_norm < 1e-20)
		return;

	Vector3d e = this->LD.F.normalized();
	double dx, F0, F1;
	F0 = 0.0;
	dx = dx0;

	// �܂��{�[���S�̂𓮂����D�i�ړ��ʂ͎w���֐������j
	for (int i = 0; i < 1e3; i++) {
		// �{�[���S�̂��׏d�����ɓ������D
		for (int i = 0; i < this->nCC; i++)
			for (int j = 0; j < this->CC[i].nBL; j++)
				this->CC[i].BL[j].x += e * dx;

		// �i�b�g�ɂ�����׏d�����߂�D
		Vector3d Fn_sum = Vector3d::Zero();
		for (int j = 0; j < this->nP; j++) {
			Vector2d Fbn;
			Vector3d Fnb, Tnb;
			this->BNP[j].get_F0(Fbn, Fnb, Tnb);
			Fn_sum += Fnb;
		}
		// �׏d���X�V�D
		F1 = Fn_sum.norm();

		// �׏d���O���׏d�𒴉߂��Ă���ꍇ�C���[�v�𔲂���D
		if (F1 > F_norm) {
			// ���`��Ԃł����悻�̈ʒu�ɖ߂��D

			for (int i = 0; i < this->nCC; i++)
				for (int j = 0; j < this->CC[i].nBL; j++)
					this->CC[i].BL[j].x -= (F1 - F_norm) / (F1 - F0) * e * dx;
			break;
		}
		// �␳�l���X�V�D
		dx *= 1.1;
		F0 = F1;
	}
	// ���ɃV���t�g�𓮂����D�i�ړ��ʂ͎w���֐������j
	F0 = 0.0;
	dx *= 0.1;
	for (int i = 0; i < 1e3; i++) {
		// �V���t�g���׏d�����ɓ������D
		this->ST.x += e * dx;
		this->ST.set_dx();

		// �V���t�g�ɂ�����׏d�����߂�D
		Vector3d Fs_sum = Vector3d::Zero();
		for (int j = 0; j < this->nP; j++) {
			Vector2d Fbs;
			Vector3d Fsb, Tsb;
			this->BSP[j].get_F0(Fbs, Fsb, Tsb);
			Fs_sum += Fsb;
		}
		// �׏d���X�V�D
		F1 = Fs_sum.norm();

		// �׏d���O���׏d�𒴉߂��Ă���ꍇ�C���[�v�𔲂���D
		if (F1 > F_norm) {
			// ���`��Ԃł����悻�̈ʒu�ɖ߂��D
			this->ST.x -= (F1 - F_norm) / (F1 - F0) * e * dx;
			this->ST.set_dx();
			break;
		}
		// �␳�l���X�V�D
		dx *= 1.1;
		F0 = F1;
	}
	return;
}

// �V���t�g�����[�����g�����ɉ�]���āC�ʂƐڐG����悤�ɏ����ʒu������
void BS_BallScrew::preset_y0_T(double dth0) {

	// �׏d���\���������ꍇ�C���̃��\�b�h�͎g���Ȃ��D
	Vector2d T2d(this->LD.T.y(), this->LD.T.z());
	double T_norm = 2 * T2d.norm();
	if (T_norm < 1e-20)
		return;

	Vector2d e = T2d.normalized();
	double dth, T0, T1;
	T0 = 0.0;
	dth = dth0;

	// �V���t�g���ړ��i�ړ��ʂ͎w���֐��I�ɑ����j
	for (int i = 0; i < 1e3; i++) {

		Vector3d ax = this->ST.get_ax();
		this->ST.set_ax(Vector3d(1.0, ax.y() + dth * e[1], ax.z() - dth * e[0]));
		this->ST.set_dx();

		Vector2d Ts_sum = Vector2d::Zero();
		for (int j = 0; j < this->nP; j++) {
			Vector2d Fbs;
			Vector3d Fsb, Tsb;
			this->BSP[j].get_F0(Fbs, Fsb, Tsb);
			Ts_sum += Vector2d(Tsb.y(), Tsb.z());
		}
		// �׏d���X�V�D
		T1 = Ts_sum.norm();

		// �׏d���O���׏d�𒴉߂��Ă���ꍇ�C���[�v�𔲂���D
		if (T1 > T_norm) {
			// ���`��Ԃł����悻�̈ʒu�ɖ߂��D
			double dth_ = -(T1 - T_norm) / (T1 - T0) * dth;
			Vector3d ax = this->ST.get_ax();
			this->ST.set_ax(Vector3d(1.0, ax.y() + dth_ * e[1], ax.z() - dth_ * e[0]));
			this->ST.set_dx();
			break;
		}
		// �␳�l���X�V�D
		dth *= 1.1;
		T0 = T1;
	}
	return;
}

// �i�b�g���Ƀ{�[�����G���܂ł������郁�\�b�h�D
void BS_BallScrew::preset_y0_x(double dx) {

	for (int i = 0; i < this->nP; i++) {
		// �܂������̂��Ƃ̈ʒu�ɖ߂��Ă��D
		Vector3d eta_zeta = this->BNP[i].get_eta();
		double eta = eta_zeta[1];
		double zeta = eta_zeta[2];
		this->BNP[i].set_eta0(Vector2d(eta, zeta));

		for (int j = 0; j < 1e3; j++) {
			Vector2d Fbn;
			Vector3d Fnb, Tnb;
			this->BNP[i].get_F0(Fbn, Fnb, Tnb);
			double F = abs(Fbn.x());

			// �ڂ��Ă����甲����D
			if (F > 1e-20)
				break;

			eta -= dx;
			dx *= 1.1;
			this->BNP[i].set_eta0(Vector2d(eta, zeta));
		}
	}
	return;
}

// step0�i�����v�Z�j�p�̃C���^�t�F�C�X�D
void BS_BallScrew::get_y0(double*y0) {

	this->ST.get_y0(y0);

	for (int i = 0; i < this->nP; i++) {
		int i2 = i * 2;
		Vector3d eta = this->BNP[i].get_eta();
		y0[i2 + 5] = eta[1] / Rigid::l;
		y0[i2 + 6] = eta[2] / Rigid::l;
	}
	return;
}

// step0�i�����v�Z�j�F�ϐ�y0�����
void BS_BallScrew::set_y0(
	const double*y0,	// in :[-]		: �ϐ�y0
	double v0,			// in :[m/s]	: �V���t�g�i�s���x
	double w0			// in :[rad/s]	: �V���t�g��]���x
) {

	this->ST.set_y0(y0, v0, w0);

	for (int i = 0; i < this->nP; i++) {
		int i2 = i * 2;
		Vector2d eta = Vector2d(
			y0[i2 + 5] * Rigid::l,
			y0[i2 + 6] * Rigid::l
		);
		this->BNP[i].set_eta0(eta);
	}
	return;
}

// �׏d��z��ɂ��ĕԂ����\�b�h�D
void BS_BallScrew::get_F0(double*F0) {

	Vector3d Fs = this->LD.F;
	Vector3d Ts = this->LD.T;

	for (int ib = 0; ib < this->nP; ib++) {

		Vector2d Fbn, Fbs;
		Vector3d Fnb, Tnb, Fsb, Tsb;

		this->BNP[ib].get_F0(Fbn, Fnb, Tnb);
		this->BSP[ib].get_F0(Fbs, Fsb, Tsb);

		Vector2d Fb = Fbn + Fbs;

		int i2 = ib * 2;

		// �ア�o�l�i�˗́j�Dx^-1 �̌`�ɂ��邱�Ƃɂ���āC�O���Ɏ������₷������D
		double x = Vector2d(this->BNP[ib].BL->x[1], this->BNP[ib].BL->x[2]).norm();
		double k = 1e-20;

		F0[i2 + 5] = Fb[0] - k / x;
		F0[i2 + 6] = Fb[1];

		Fs += Fsb;
		Ts += Tsb;
	}
	Vector3d Ts_ = Ts / Rigid::l;

	return;
}

// �׏d��z��ɂ��ĕԂ����\�b�h�D
void BS_BallScrew::get_F1(double*F1) {

	Vector3d Fs = this->LD.F;
	Vector3d Ts = this->LD.T;

	bool is_RightHand = (this->NT->w.x() - this->ST.w.x()) > 0;	// �V���t�g���猩�ċʂ������v���Ɍ��]���Ă���Ƃ�true

	for (int ib = 0; ib < this->nP; ib++) {

		Vector2d Fbn, Fbs;
		Vector3d Fnb, Fsb, Tnb, Tsb;
		this->BNP[ib].get_F1(is_RightHand, Fbn, Fnb, Tnb);
		this->BSP[ib].get_F1(!is_RightHand, Fbs, Fsb, Tsb);
		Vector2d Fb = Fbn + Fbs;
		//double Fbn_ = Fbn.norm();
		//double Fbs_ = Fbs.norm();

		//double sin_th = (Fbn[1] * Fbs[0] - Fbn[0] * Fbs[1]) / (Fbn_ * Fbs_);

		int i2 = ib * 2;
		F1[i2 + 5] = Fb[0]; // Fbn_ - Fbs_;
		F1[i2 + 6] = Fb[1]; // sin_th * sqrt(Fbn_ * Fbs_);

		Fs += Fsb;
		Ts += Tsb;
	}
	Vector3d Ts_ = Ts / Rigid::l;

	F1[0] = Fs[0];
	F1[1] = Fs[1];
	F1[2] = Fs[2];
	F1[3] = Ts_[1];
	F1[4] = Ts_[2];

	return;
}

// �e�]���̂̏��]���葬�x�����߂郁�\�b�h�D
void BS_BallScrew::pure_Rolling(void) {

	// �S�Ă̓]���̂����]�����Ԃɐݒ肷��D
	for (int i = 0; i < this->nP; i++) {

		// �e�a�ł̐ڐG�_�������߂�D
		double dxn[2], dxs[2]; Vector3d pn[2], ps[2], v_guide;
		int nn = this->BNP[i].num_Contact(dxn, pn);
		int ns = this->BSP[i].num_Contact(dxs, ps);

		// �e���ڋߗʂ��傫��������׏d����D
		Vector3d pn0 = (dxn[0] > dxn[1]) ? pn[0] : pn[1];
		Vector3d ps0 = (dxs[0] > dxs[1]) ? ps[0] : ps[1];

		// �ڐG�_�������������ʐi�s�����Ƃ���D�����ꍇ�͒��Ԃ����D
		Vector3d ev0 =
			(nn > ns) ? this->BNP[i].get_e() :
			(nn < ns) ? this->BSP[i].get_e() :
			(this->BNP[i].get_e() + this->BSP[i].get_e()).normalized();

		// �ʐi�s���������Ō����C�e��׏d�\�ʂ̕����t�����x���Z�o����D
		double vn0 = this->BNP[i].CY->surface_velocity(pn0).dot(ev0);
		double vs0 = this->BSP[i].CY->surface_velocity(ps0).dot(ev0);

		double r0 = (pn0 - ps0).norm() / 2;
		Vector3d ep0 = (pn0 - ps0) / (2 * r0);
		Vector3d ew0 = ev0.cross(ep0).normalized();
		double v0 = (vn0 + vs0) / 2;
		double w0 = (v0 - vn0) / r0;

		this->BNP[i].BL->v = v0 * ev0;
		this->BNP[i].BL->w = w0 * ew0;

		this->BNP[i].stt.ev0 = ev0;
		this->BNP[i].stt.w0 = w0 * ew0;
	}
	return;
}

bool BS_BallScrew::get_y1(int ib, double * y1) {
	
	Vector3d yx = this->BNP[ib].get_eta() / Rigid::l;
	double   yv = this->BNP[ib].stt_get_v0() * Rigid::t / Rigid::l;
	Vector3d yw = this->BNP[ib].stt_get_dw0() * Rigid::t;

	y1[0] = yx[1];
	y1[1] = yx[2];
	y1[2] = yv;
	y1[3] = yw[0];
	y1[4] = yw[1];
	y1[5] = yw[2];

	bool has_next = ib < (this->nP - 1);

	return has_next;
}

void BS_BallScrew::set_y1(int ib, const double * y1) {

	Vector2d x = Vector2d(y1[0], y1[1]) * Rigid::l;
	double   v = y1[2] / Rigid::t * Rigid::l;
	Vector3d w = Vector3d(y1[3], y1[4], y1[5])
		/ Rigid::t;

	this->BNP[ib].set_eta0(x);
	this->BNP[ib].stt_set_v0(v);
	this->BNP[ib].stt_set_dw0(w);

	return;
}

void BS_BallScrew::get_F1(int ib, double * f1) {

	Vector3d Fbn, Tbn, Tnb, Fnbs, Tnbs;
	this->BNP[ib].get_FT(Fbn, Tbn, Tnb, Fnbs, Tnbs);

	Vector3d Fbs, Tbs, Tsb, Fsbs, Tsbs;
	this->BSP[ib].get_FT(Fbs, Tbs, Tsb, Fsbs, Tsbs);

	Vector3d BL_F = Fbn + Fbs + this->BNP[ib].BL->get_mg();
	Vector3d BL_T = Tbn + Tbs;

	for (int j = 0; j < 3; j++) {
		f1[j + 0] = BL_F[j];
		f1[j + 3] = BL_T[j];
	}
	return;
}

void BS_BallScrew::get_y2(double * y2) {

	this->ST.get_y0(y2);

	for (int ib = 0; ib < this->nP; ib++) {

		int i6 = ib * 6;
		Vector3d yx = this->BNP[ib].get_eta() / Rigid::l;
		double   yv = this->BNP[ib].stt_get_v0() * Rigid::t / Rigid::l;
		Vector3d yw = this->BNP[ib].stt_get_dw0() * Rigid::t;

		y2[i6 + 5] = yx[1];
		y2[i6 + 6] = yx[2];
		y2[i6 + 7] = yv;
		y2[i6 + 8] = yw[0];
		y2[i6 + 9] = yw[1];
		y2[i6 + 10] = yw[2];
	}
	return;
}


void BS_BallScrew::set_y2(const double * y2, double v0, double w0) {

	this->ST.set_y0(y2, v0, w0);

	for (int ib = 0; ib < this->nP; ib++) {

		int i6 = ib * 6;
		Vector2d x = Vector2d(y2[i6 + 5], y2[i6 + 6]) * Rigid::l;
		double   v = y2[i6 + 7] / Rigid::t * Rigid::l;
		Vector3d w = Vector3d(y2[i6 + 8], y2[i6 + 9], y2[i6 + 10])
			/ Rigid::t;

		this->BNP[ib].set_eta0(x);
		this->BNP[ib].stt_set_v0(v);
		this->BNP[ib].stt_set_dw0(w);
	}
	return;
}

void BS_BallScrew::get_F2(double * f2) {

	// �e�ʂ̐ڐG�͂����߂�D
	Vector3d ST_F = this->LD.F + this->ST.get_mg();
	Vector3d ST_T = this->LD.T;

	for (int i = 0; i < this->nP; i++) {

		Vector3d Fbn, Tbn, Tnb, Fnbs, Tnbs;
		this->BNP[i].get_FT(Fbn, Tbn, Tnb, Fnbs, Tnbs);

		Vector3d Fbs, Tbs, Tsb, Fsbs, Tsbs;
		this->BSP[i].get_FT(Fbs, Tbs, Tsb, Fsbs, Tsbs);
		ST_F += -Fbs;	// �ʂɂ�����͂̔���p�ɂ�蕄�����]�D
		ST_T += Tsb;

		Vector3d BL_F = Fbn + Fbs + this->BNP[i].BL->get_mg();
		Vector3d BL_T = Tbn + Tbs;

		for (int j = 0; j < 3; j++) {
			f2[5 + 6 * i + j] = BL_F[j];
			f2[8 + 6 * i + j] = BL_T[j];
		}
	}
	double y2[5];
	this->ST.get_y0(y2);

	for (int i = 0; i < 3; i++)
		f2[i] = (this->ST.x_const[i]) ? y2[i] : ST_F[i];

	f2[3] = (this->ST.Rx_const[1]) ? y2[3] : ST_T[1];
	f2[4] = (this->ST.Rx_const[2]) ? y2[4] : ST_T[2];

	return;
}


void BS_BallScrew::init_dyn0(void) {

	this->ST.init_dyn0();

	for (int ib = 0; ib < this->nP; ib++) {
		this->BNP[ib].mem_BLv = this->BNP[ib].BL->v;
		this->BNP[ib].BL->v.setZero();
	}
	return;
}

void BS_BallScrew::deinit_dyn0(double v0, double w0) {

	this->ST.deinit_dyn0(v0, w0);

	for (int ib = 0; ib < this->nP; ib++)
		this->BNP[ib].BL->v = this->BNP[ib].mem_BLv;

	return;
}

void BS_BallScrew::get_dyn_y0(double * y0) {

	this->ST.get_dyn_y0(y0);

	for (int ib = 0; ib < this->nP; ib++) {
		int i4 = ib * 4;
		Vector3d eta = this->BNP[ib].get_eta0()
			/ Rigid::l;
		y0[i4 + 10] = eta[1];
		y0[i4 + 11] = eta[2];

		Vector3d etav = this->BNP[ib].get_etav0()
			/ Rigid::l * Rigid::t;
		y0[i4 + 12] = etav[1];
		y0[i4 + 13] = etav[2];
	}
	return;
}

// step0�i�����v�Z�j�F�ϐ�y0�����
void BS_BallScrew::set_dyn_y0(const double*y0) {

	this->ST.set_dyn_y0(y0);

	for (int ib = 0; ib < this->nP; ib++) {
		int i4 = ib * 4;

		Vector2d eta = Vector2d(
			y0[i4 + 10], y0[i4 + 11]
		)* Rigid::l;

		this->BNP[ib].set_eta0(eta);

		Vector2d etav = Vector2d(
			y0[i4 + 12], y0[i4 + 13]
		)* Rigid::l / Rigid::t;

		this->BNP[ib].set_etav0(etav);
	}
	return;
}

// �׏d��z��ɂ��ĕԂ����\�b�h�D
void BS_BallScrew::get_dyn_dydt0(double*dydt) {

	Vector3d Fst = this->LD.F;
	Vector3d Tst = this->LD.T;

	bool is_RightHand = (this->NT->w.x() - this->ST.w.x()) > 0;	// �V���t�g���猩�ċʂ������v���Ɍ��]���Ă���Ƃ�true

	for (int ib = 0; ib < this->nP; ib++) {

		int i4 = ib * 4;

		Vector3d etavn, etavs, Fbn, Fbs, Fnb, Fsb, Tnb, Tsb;
		this->BNP[ib].get_dyn_F0(is_RightHand, etavn, Fbn, Fnb, Tnb);
		this->BSP[ib].get_dyn_F0(!is_RightHand, etavs, Fbs, Fsb, Tsb);

		//std::cout << "bef" << etavn << std::endl << std::endl;

		etavn = this->BNP[ib].get_etav0();

		//std::cout << "aft" << etavn << std::endl << std::endl;

		dydt[i4 + 10] = etavn[1] / Rigid::l * Rigid::t;
		dydt[i4 + 11] = etavn[2] / Rigid::l * Rigid::t;

		Matrix3d xyz2eta = this->BNP[ib].get_xyz2eta();
		Vector3d aeta = this->BNP[ib].CY->to_etavector(
			Fbn + Fbs, xyz2eta)
			* this->BNP[ib].BL->m_inv
			/ Rigid::l * Rigid::t * Rigid::t;

		dydt[i4 + 12] = aeta[1];
		dydt[i4 + 13] = aeta[2];

		Fst += Fsb;
		Tst += Tsb;
	}
	Vector3d    dxdt = this->ST.get_dxdt();
	Vector3d    dvdt = this->ST.get_dvdt(Fst);
	Quaterniond dqdt = this->ST.get_dqdt();
	Vector3d    dwdt = this->ST.get_dwdt(Tst);

	for (int i = 0; i < 3; i++) {
		dydt[i + 0] = dxdt[i] / Rigid::l * Rigid::t;
		dydt[i + 3] = dvdt[i] / Rigid::l * Rigid::t * Rigid::t;
	}
	dydt[6] = dqdt.y() * Rigid::t;
	dydt[7] = dqdt.z() * Rigid::t;
	dydt[8] = dwdt.y() * Rigid::t * Rigid::t;
	dydt[9] = dwdt.z() * Rigid::t * Rigid::t;

	//for (size_t i = 90; i < 94; i++)
	//	std::cout << "\t" << dydt[i];

	//std::cout << std::endl << std::endl;
	return;
}

void BS_BallScrew::set_dyn_y1(const double*y) {

	Vector3d x, v, w; Quaterniond q;
	x.x() = y[0];
	x.y() = y[1];
	x.z() = y[2];
	v.x() = y[3];
	v.y() = y[4];
	v.z() = y[5];
	q.w() = y[6];
	q.x() = y[7];
	q.y() = y[8];
	q.z() = y[9];
	w.x() = y[10];
	w.y() = y[11];
	w.z() = y[12];
	this->NT->set_y_(x, v, q, w);

	x.x() = y[13];
	x.y() = y[14];
	x.z() = y[15];
	v.x() = y[16];
	v.y() = y[17];
	v.z() = y[18];
	q.w() = y[19];
	q.x() = y[20];
	q.y() = y[21];
	q.z() = y[22];
	w.x() = y[23];
	w.y() = y[24];
	w.z() = y[25];
	this->ST.set_y_(x, v, q, w);

	for (int i = 0; i < this->nP; i++) {
		int i13 = i * 13 + 26;
		x.x() = y[i13 + 0];
		x.y() = y[i13 + 1];
		x.z() = y[i13 + 2];
		v.x() = y[i13 + 3];
		v.y() = y[i13 + 4];
		v.z() = y[i13 + 5];
		q.w() = y[i13 + 6];
		q.x() = y[i13 + 7];
		q.y() = y[i13 + 8];
		q.z() = y[i13 + 9];
		w.x() = y[i13 + 10];
		w.y() = y[i13 + 11];
		w.z() = y[i13 + 12];
		this->BNP[i].BL->set_y(x, v, q, w);
	}
	return;
}

void BS_BallScrew::get_dyn_y1(double*y) {

	Vector3d x, v, w; Quaterniond q;
	this->NT->get_y(x, v, q, w);
	y[0] = x.x();
	y[1] = x.y();
	y[2] = x.z();
	y[3] = v.x();
	y[4] = v.y();
	y[5] = v.z();
	y[6] = q.w();
	y[7] = q.x();
	y[8] = q.y();
	y[9] = q.z();
	y[10] = w.x();
	y[11] = w.y();
	y[12] = w.z();

	this->ST.get_y(x, v, q, w);
	y[13] = x.x();
	y[14] = x.y();
	y[15] = x.z();
	y[16] = v.x();
	y[17] = v.y();
	y[18] = v.z();
	y[19] = q.w();
	y[20] = q.x();
	y[21] = q.y();
	y[22] = q.z();
	y[23] = w.x();
	y[24] = w.y();
	y[25] = w.z();

	for (int i = 0; i < this->nP; i++) {
		int i13 = i * 13 + 26;
		this->BNP[i].BL->get_y(x, v, q, w);
		y[i13 + 0] = x.x();
		y[i13 + 1] = x.y();
		y[i13 + 2] = x.z();
		y[i13 + 3] = v.x();
		y[i13 + 4] = v.y();
		y[i13 + 5] = v.z();
		y[i13 + 6] = q.w();
		y[i13 + 7] = q.x();
		y[i13 + 8] = q.y();
		y[i13 + 9] = q.z();
		y[i13 + 10] = w.x();
		y[i13 + 11] = w.y();
		y[i13 + 12] = w.z();
	}
	return;
}

void BS_BallScrew::get_dyn_dydt1(
	double*dydt,	// out:	�����z��
	double dvdt,	// in:	
	double dwdt 	// in:	
) {
	// �e�ʂ̐ڐG�͂����߂�D
	this->NT->set_Fall_Zero();
	this->ST.set_Fall_Zero();

	for (int i = 0; i < this->nP; i++) {

		Vector3d Fbn, Tbn, Tnb, Fnbs, Tnbs;
		this->BNP[i].get_FT(Fbn, Tbn, Tnb, Fnbs, Tnbs);
		this->BNP[i].CY->F += -Fbn;
		this->BNP[i].CY->T += Tnb;
		this->BNP[i].CY->Fs += Fnbs;
		this->BNP[i].CY->Ts += Tnbs;

		Vector3d Fbs, Tbs, Tsb, Fsbs, Tsbs;
		this->BSP[i].get_FT(Fbs, Tbs, Tsb, Fsbs, Tsbs);
		this->BSP[i].CY->F += -Fbs;
		this->BSP[i].CY->T += Tsb;
		this->BSP[i].CY->Fs += Fsbs;
		this->BSP[i].CY->Ts += Tsbs;

		// �d�͉����x�����Z�D
		this->BNP[i].BL->F = Fbn + Fbs + this->BNP[i].BL->get_mg();
		this->BNP[i].BL->T = Tbn + Tbs;
		double dydt_[13];
		this->BNP[i].BL->get_dydt(this->BNP[i].BL->F, this->BNP[i].BL->T, dydt_);

		// ���[�v���Ƃɖ߂�l�ɑ���D
		int i13 = i * 13 + 26;
		for (int j = 0; j < 13; j++)
			dydt[i13 + j] = dydt_[j];

		this->NT->F += -Fbn;
		this->NT->T += Tnb;
		this->ST.F += -Fbs;
		this->ST.T += Tsb;
	}
	// �ꉞ�������Ƃ邪�CNT��v_const=true�Ȃ̂�0���߂�l�ƂȂ�D
	double dydtn[13];
	this->NT->get_dydt(this->NT->F, this->NT->T, dydtn);
	for (int j = 0; j < 13; j++)
		dydt[j] = dydtn[j];

	Vector3d ST_F;
	ST_F = this->ST.F + this->LD.F + this->ST.get_mg();

	Vector3d ST_T;
	ST_T = this->ST.T + this->LD.T;

	double dydts[13];
	this->ST.get_dydt_(ST_F, ST_T, dvdt, dwdt, dydts);

	for (int j = 0; j < 13; j++)
		dydt[j + 13] = dydts[j];

	return;
}

//void BS_BallScrew::Shaft_Lock(void) {
//
//	this->ST.set_const(true, true);
//
//	return;
//};

// �O���׏d����
void BS_BallScrew::set_load(double *F, double *T) {
	for (int i = 0; i < 3; i++) {
		this->LD.F[i] = F[i];
		this->LD.T[i] = T[i];
	}
	return;
};

// ���݂̏�Ԃ��O���ϐ��ɓ��́D
void BS_BallScrew::save(BS_Out&OUT) {

	for (int i = 0; i < this->nCC; i++)
		for (int j = 0; j < this->CC[i].nBL; j++)
			this->CC[i].BL[j].save(OUT.CC[i].BL[j].x, OUT.CC[i].BL[j].v, OUT.CC[i].BL[j].q, OUT.CC[i].BL[j].w, OUT.CC[i].BL[j].ax, OUT.CC[i].BL[j].F, OUT.CC[i].BL[j].T);

	this->NT->save(OUT.NT.x, OUT.NT.v, OUT.NT.q, OUT.NT.w, OUT.NT.ax, OUT.NT.F, OUT.NT.T);

	for (int i = 0; i < this->NT->nCY; i++)
		this->NT->CY[i].save(OUT.NT_CY[i]);

	this->ST.save(OUT.ST.x, OUT.ST.v, OUT.ST.q, OUT.ST.w, OUT.ST.ax, OUT.ST.F, OUT.ST.T);

	for (int i = 0; i < this->ST.nCY; i++)
		this->ST.CY[i].save(OUT.ST_CY[i]);

	for (int i = 0; i < this->nP; i++)
		this->BNP[i].save(OUT.BNP[i]);

	for (int i = 0; i < this->nP; i++)
		this->BSP[i].save(OUT.BSP[i]);

	return;
}

BS_BallScrew::BS_BallScrew() {
	this->CC = NULL;
	this->BBP = NULL;
	this->BNP = NULL;
	this->BSP = NULL;
	return;
}

BS_BallScrew::~BS_BallScrew() {

	if (this->CC != NULL)
		delete[] this->CC;
	if (this->BBP != NULL)
		delete[] this->BBP;
	if (this->BNP != NULL)
		delete[] this->BNP;
	if (this->BSP != NULL)
		delete[] this->BSP;
	return;
}













//// ���S�͂ƃW���C�����[�����g�D�Ƃ肠�����x�^�����ŁD
//Vector3d x = this->BL[ib].x;
//Vector3d v = this->BL[ib].v;
//Vector3d w = this->BL[ib].w;
//double rv = x[1] * v[2] - x[2] * v[1];
//Vector3d x_ = Vector3d(0.0, x[1], x[2]);
//Vector3d e = x_.normalized();
//double r = x_.norm();
//double Fc = this->BL[ib].m * rv * rv / r / r / r;
//Vector3d W = Vector3d(rv / r / r, 0.0, 0.0);
//Vector3d L = this->BL[ib].I.cwiseProduct(w);
//Vector3d Tc = W.cross(L);

//// stepp�i�^���v�Z�j�p�̃C���^�t�F�C�X�D
//void BS_BallScrew::get_yp(double*yp) {
//
//	for (int i = 0; i < 2; i++)
//		yp[i] = this->NT->CY[i].x.x();
//
//	for (int i = 0; i < this->nBL; i++) {
//		int i2 = i * 2;
//		Vector3d eta = this->BSP[i].get_eta();
//		yp[i2 + 2] = eta[1];
//		yp[i2 + 3] = eta[2];
//	}
//	return;
//}
//
//// stepp�i�^���v�Z�j�p�̃C���^�t�F�C�X�D
//void BS_BallScrew::set_yp(const double*yp) {
//
//	for (int i = 0; i < 2; i++)
//		this->NT->dx[i] = Vector3d(yp[i], 0.0, 0.0);
//	this->NT->set_dx();
//
//	for (int i = 0; i < this->nBL; i++) {
//		int i2 = i * 2;
//		Vector2d eta = Vector2d(yp[i2 + 2], yp[i2 + 3]);
//		this->BSP[i].set_etap(eta);
//	}
//	return;
//}
//
//// �׏d��z��ɂ��ĕԂ����\�b�h�D
//void BS_BallScrew::get_Fp(double*Fp) {
//
//	this->NT->CY[0].F = Vector3d(-this->F_pre, 0.0, 0.0);
//	this->NT->CY[1].F = Vector3d(this->F_pre, 0.0, 0.0);
//
//	for (int i = 0; i < this->nBL; i++) {
//
//		Vector2d Fbn, Fbs;
//		Vector3d Fnb, Tnb, Fsb, Tsb;
//
//		this->BNP[i].get_F0(Fbn, Fnb, Tnb);
//		this->BSP[i].get_F0(Fbs, Fsb, Tsb);
//
//		Vector2d Fb = Fbn + Fbs;
//
//		int i2 = i * 2;
//		Fp[i2 + 2] = Fb[0];
//		Fp[i2 + 3] = Fb[1];
//
//		this->BNP[i].CY->F += Fnb;
//	}
//	for (int i = 0; i < 2; i++)
//		Fp[i] = this->NT->CY[i].F.x();
//
//	return;
//}

//
//for (size_t i = 0; i < this->nBL - 1; i++) {
//	Vector3d hoge = this->BNP[i].get_e();
//	Vector3d fuga = this->BL[i+1].x - this->BL[i].x;
//	cout << hoge.dot(fuga) / fuga.norm() << endl;
//}
//
	//for (int i = 0; i < this->nBL; i++) {
	//	Vector3d e = this->BNP[i].get_e();
	//	this->BNP[i].set_e(e);
	//	this->BSP[i].set_e(e);
	//}


//// step1�i�����v�Z�j�p�̃C���^�t�F�C�X�D
//void BS_BallScrew::get_y1(double*y1) {
//
//	int ib = this->ib;
//
//	Vector3d eta = this->BNP[ib].get_eta();
//	y1[0] = eta[1] / Rigid::l;
//	y1[1] = eta[2] / Rigid::l;
//
//	return;
//}
//
//void BS_BallScrew::set_y1(const double*y1) {
//
//	int ib = this->ib;
//
//	Vector2d eta = Vector2d(
//		y1[0] * Rigid::l,
//		y1[1] * Rigid::l
//	);
//	this->BNP[ib].set_eta0(eta);
//
//	return;
//}


//// step2�i�����E���C�v�Z�j�p�̃C���^�t�F�C�X�D
//void BS_BallScrew::set_y2(const double*y2) {
//
//	int ib = this->ib;
//
//	Vector2d eta = Vector2d(
//		y2[0] * Rigid::l,
//		y2[1] * Rigid::l
//	);
//	this->BNP[ib].set_eta0(eta);
//
//	Vector3d v = Vector3d(y2[2], y2[3], y2[4]) / Rigid::t * Rigid::l;
//
//	this->BL[ib].v = v;
//
//	Vector3d w = Vector3d(y2[5], y2[6], y2[7]) / Rigid::t;
//	this->BL[ib].w = w;
//
//	return;
//}
//
//// �׏d��z��ɂ��ĕԂ����\�b�h�D
//void BS_BallScrew::get_F2(double*F2) {
//
//	int ib = this->ib;
//
//	Vector3d Fn, Tn, Fs, Ts, T_n, T_s;
//	this->BNP[ib].get_FT(Fn, Tn, T_n);
//	this->BSP[ib].get_FT(Fs, Ts, T_s);
//
//	// ���S�͂ƃW���C�����[�����g�D�Ƃ肠�����x�^�����ŁD
//	Vector3d x = this->BL[ib].x;
//	Vector3d v = this->BL[ib].v;
//	Vector3d w = this->BL[ib].w;
//	double rv = x[1] * v[2] - x[2] * v[1];
//	Vector3d x_ = Vector3d(0.0, x[1], x[2]);
//	Vector3d e = x_.normalized();
//	double r = x_.norm();
//	double Fc = this->BL[ib].m * rv * rv / r / r / r;
//	Vector3d W = Vector3d(rv / r / r, 0.0, 0.0);
//	Vector3d L = this->BL[ib].I.cwiseProduct(w);
//	Vector3d Tc = W.cross(L);
//
//	Vector3d Fb = Fn + Fs + e * Fc;
//	Vector3d Tb = (Tn + Ts + Tc) / Rigid::l;
//
//	F2[0] = Fb[0];
//	F2[1] = Fb[1];
//	F2[2] = Fb[2];
//	F2[3] = Tb[0];
//	F2[4] = Tb[1];
//	F2[5] = Tb[2];
//
//	F2[6] = 0.0;
//	F2[7] = 0.0;
//
//	return;
//}

//// step2�i�����E���C�v�Z�j�p�̃C���^�t�F�C�X�D
//void BS_BallScrew::get_y2(double*y2) {
//
//	int ib = this->ib;
//
//	Vector3d eta = this->BNP[ib].get_eta();
//	y2[0] = eta[1] / Rigid::l;
//	y2[1] = eta[2] / Rigid::l;
//
//	//y2[0] = 0.011e-3;
//	//y2[1] = 0.020e-3;
//
//	Vector3d v = this->BL[ib].v * Rigid::t / Rigid::l;
//	y2[2] = v[0];
//	y2[3] = v[1];
//	y2[4] = v[2];
//
//	Vector3d w = this->BL[ib].w * Rigid::t;
//	y2[5] = w[0];
//	y2[6] = w[1];
//	y2[7] = w[2];
//
//	return;
//}


		//std::cout
		//	<< i << std::endl
		//	<< nn << std::endl
		//	<< ns << std::endl
		//	<< this->BNP[i].get_us(pn0) << std::endl
		//	<< this->BSP[i].get_us(ps0) << std::endl
		//	<< std::endl;
