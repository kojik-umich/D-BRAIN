#pragma once

class B4P_In {
public:
	virtual void Nothing(void)=0;	// ���ۃN���X�ł��邱�Ƃ��������߁C�Ӗ��̂Ȃ��֐����`����D
	
public:
	int ballnum;					// �]���̂̌�
	double balldia;				// ��{�ʌa�i�����l�v�Z�݂̂ɗp����j
	double ballpcd;				// �]����pcd�i�w4�̍aR���Spcd�̕��ϒl�x�����D�����l�v�Z�݂̂Ɏg�p�j
	double cos_alp0;			// �ڐG�p�̗]���i�����l�v�Z�݂̂Ɏg�p�j

	struct Tribology {
		enum RollingResistance {
			RollingResistanceNothing = 0, // �]����S����R�Ȃ�
			GoksemAihara = 1,	// Goksem-�����̎��ibaltac�������j
			Aihara = 2,			// �����̎�
			Fujiwara = 3,		// �����̎�
			Houpert = 4,		// Houpert�̎�
		} rollingresistance;

		enum Coulomb {
			CoulombNothing = 0,	// �N�[�������C�Ȃ��D
			Tangent = 1,	// �^���W�F���g�J�[�u�D
		} coulomb;
		double coulomb_slope;

		enum FilmThickness {
			FilmThicknessNothing = 0,	// �����Ȃ��D
			HamrockDowsonHc = 1,			// HamrockDowson���̒������������𖀎C�v�Z�Ɏg�p�D
			HamrockDowsonHmin = 2,			// HamrockDowson���̍ŏ����������𖀎C�v�Z�Ɏg�p�D
		} filmThickness;

		enum Hysteresis {
			HysteresisNothing = 0,		// �q�X�e���V�X�Ȃ��D
			Kakuta = 1,					// �p�c�̎�
		}hysteresis;
		double hysteresis_factor;		// �q�X�e���V�X�����W��
	} TB;

	enum CageType {
		spherical_machined_cage		= 0,	// ���ʂ��ݔ���
		snap_cage					= 1,	// ���^
		cylindrical_machined_cage	= 2,	// �~�����ݔ���
		corrugated_press_cage		= 3,	// �g�^�v���X
		cylindrical_snap_cage		= 4,	// �~�����^TY
		snap_press_cage				= 5,	// ���^�v���X
		full_ball					= -2,	// ���ʎ���
		no_cage						= -1	// �ێ��햳��
	} cage_type;
	   
	// �]���̖̂��x�C�����C�ʌa�Ȃǂ͋ʔԍ����Ƃɒ�`
	struct Ball {
		double den;				// �]���̖��x[kg/m^3]
		double E;				// �]���̍���[Pa]
		double por;				// �]���̃|�A�\����
		double rms;				// �]���̑e��rms
		double dia;				// �]���̊O�a�i�]���̖��ɒ�`�j
	} *BL;
									   

									   
	struct Ring {
		double den;				// ���O�֖��x[kg/m^3]
		double E;				// ���O�փ����O��[Pa]
		double por;				// ���O�փ|�A�\����
		double rms;				// ���O�֑e��rms
		double R[2];			// ���O�֍aR���a[m] 0:-x���a, 1:+x���a
		double Rox[2];			// ���O�֍aR���Sx���W[m]
		double Rod[2];			// ���O�֍aR���S���a[m]
		double hedge[2];		// ���O�֍a������[m]
		double m;				// ���O�֏d��[kg]
		double Ix;				// ���O�֊������[�����gx����[kg�Em^2]
		double Iyz;				// ���O�֊������[�����gy����[kg�Em^2]
	} IR, OR;
	
	struct BallRingPair {
		double mu;				// ���C�W��	�i�]����-���֊ԁj
		double dzeta;			// ������		�i�]����-���֊ԁj
	} BIP, BOP;

	struct BallCagePair {
		double mu;				// ���C�W��	�i�]����-�ێ���ԁj
		double dzeta;			// ������		�i�]����-�ێ���ԁj
	} BCP;

	struct CageRingPair {
		double mu;				// ���C�W��	�i�ێ���-���֊ԁj
		double dzeta;			// ������		�i�ێ���-���֊ԁj
	}CIP, COP;

	// ���^�ێ���
	struct SnapCage {
		double dout;			// ���^�ێ���O�a[m]			// init�Ŏg��
		double din;				// ���^�ێ�����a[m]			// init�Ŏg��
		double h;				// ���^�ێ��퍂��[m]			// init�Ŏg��
		double jc;				// ���^�ێ���jc��[m] ���i�ێ����PCD�j�|�i�����PCD�j	// init�Ŏg��
		double R;				// ���^�ێ���|�P�b�g���a����[m](pocket clearance)					// init�Ŏg��
		double ropen;			// ���^�ێ���|�P�b�g�J�������a[m](pocket opening Radius)			// init�Ŏg��
		double Kface;			// �]����-�ێ��퍄���i�ʐڐG�j
		double Kedgein;			// �]����-�ێ��퍄���i�G�b�W�����j
		double Kedgeout;		// �]����-�ێ��퍄���i�G�b�W�O���j
		double Kcornerin;		// �]����-�ێ��퍄���i�p�����j
		double Kcornerout;		// �]����-�ێ��퍄���i�p�O���j
		double Kopen;			// �]����-�ێ��퍄���i�J�����j
	} Snap;

	struct Cage {
		double den;				// ���^�ێ��햧�x[kg/m^3]	// �ڐG���������v�Z�ȊO�Ŏg��Ȃ�
		double E;				// ���^�ێ��탄���O��[Pa]	// �ڐG���������v�Z�ȊO�Ŏg��Ȃ�
		double por;				// ���^�ێ���|�A�\����		// �ڐG���������v�Z�ȊO�Ŏg��Ȃ�
		double rms;				// ���^�ێ���e��rms			// ���C�W�������v�Z�ȊO�Ŏg��Ȃ�
		double rmg0[3];			// ���^�ێ���@�d�S����ɂ����􉽒��S�̈ʒu�x�N�g���iY�EZ������0�ɂ������́j
		double m;				// ���^�ێ���d��[kg]
		double Ix;				// ���^�ێ��튵�����[�����gx����[kg�Em^2]
		double Iyz;				// ���^�ێ��튵�����[�����gy�����Ez����[kg�Em^2]
	} Cage;
			
	struct Lubrication {
		double eta0;			// �S�x [Pa*s]
		double beta0;			// ���x�S�x�W�� [K^-1]
		double k0;				// ���M�`���� [W/(K*m)]
		double alpha0;			// ���͔S�x�W�� [Pa^-1]
		double lm0;				// ���j�X�J�X���� [m]
	} LB;

	double LoadIn[6];			// ���։׏dx,y,z,Rx,Ry,Rz[N][N�Em]
		   
	double omegair;				// ���։�]���x[rad/s]
	double omegaor;				// �O�։�]���x[rad/s]
	int msmax = 21;				// �ڐG�ȉ~������
	
	struct Rigid {
		double l;
		double t;
		double g[3];			// �d��(x,y,z)
	} rigid;

	struct Bound {
		bool v_const[3];
		bool w_const[3];
	} bound;
};





































	// bool autocalc_Inner_m;		// ���֏d�ʎ����v�Z�itrue:�����v�Z�Cfalse:�����v�Z���Ȃ��j
	// bool autocalc_Inner_Ix;		// ���֊������[�����gx���������v�Z�itrue:�����v�Z�Cfalse:�����v�Z���Ȃ��j
	// bool autocalc_Inner_Iy;		// ���֊������[�����gy���������v�Z�itrue:�����v�Z�Cfalse:�����v�Z���Ȃ��j
	// bool gravity_pulls_ring;	// �O���ււ̏d�͍l��(true:�d�͍l������Cfalse:�d�͍l�����Ȃ��C�f�t�H���g�Ftrue)
	// bool autocalc_Snap_m;		// ���^�ێ���@�d�ʎ����v�Z�itrue:�����v�Z�Cfalse:�����v�Z���Ȃ��j
	// bool autocalc_Snap_Ix;		// ���^�ێ���@�������[�����gx���������v�Z�itrue:�����v�Z�Cfalse:�����v�Z���Ȃ��j
	// bool autocalc_Snap_IyIz;	// ���^�ێ���@�������[�����gy�����Ez���������v�Z�itrue:�����v�Z�Cfalse:�����v�Z���Ȃ��j
	// Vector3d rmg_;			// �ێ���̏d�S���猩���􉽒��S��y�����ʒu[m]?
	// Vector3d rmg0_;			// �ێ���̏d�S���猩���􉽒��S��y�����ʒu[m]
	// bool autocalc_Snap_rmg[3];	// ���^�ێ���@�􉽒��S�����v�Z�itrue:�����v�Z�Cfalse:�����v�Z���Ȃ��j
	// DParam Outer;
	// bool autocalc_Outer_m;		// �O�֏d�ʎ����v�Z�itrue:�����v�Z�Cfalse:�����v�Z���Ȃ��j
	// bool autocalc_Outer_Ix;		// �O�֊������[�����gx���������v�Z�itrue:�����v�Z�Cfalse:�����v�Z���Ȃ��j
	// bool autocalc_Outer_Iy;		// �O�֊������[�����gy���������v�Z�itrue:�����v�Z�Cfalse:�����v�Z���Ȃ��j
	// DParam ContaBI;
	// DParam ContaBO;
	// DParam ContaBC;
	// double BC_Kinput;		// �������͒l	�i�]����-�ێ���ԁj
	// DParam ContaCI;
	// DParam ContaCO;	





