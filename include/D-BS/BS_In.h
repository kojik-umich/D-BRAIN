#pragma once
#include <vector>

class BS_In {
public:
	virtual void Nothing(void) = 0;	// ���ۃN���X�ł��邱�Ƃ��������߁C�Ӗ��̂Ȃ��֐����`����D

public:
	int ballnum;
	struct Tribology {
		enum RollingResistance {
			RollingResistanceNothing = 0,	// �]���薀�C�Ȃ�
			GoksemAihara = 1,	// �����̎�
			Aihara = 2,	// �����̎�
			Fijiwara = 3,	// �����̎�
			Houpert = 4	// Houpert�̎�
		} rollingresistance;

		enum Coulomb {
			CoulombNothing = 0,		// �N�[�������C�Ȃ��D
			Tangent = 1,			// �^���W�F���g�J�[�u�D
			SimpleCoulomb = 2
		} coulomb;
		double coulomb_slope;

		enum Ellipse {
			Mesh1d = 0,	// 1�������b�V��
			Mesh2d = 1	// 2�������b�V��
		} ellipse;
		int ellipse_mesh[2];

		enum FilmThickness {
			FilmThicknessNothing = 0,	// �����Ȃ��D
			HamrockDowsonHc = 1,			// HamrockDowson���ɂ�钆�����������D
			HamrockDowsonHmin = 2,			// HamrockDowson���ɂ��ŏ����������D
		} filmThickness;

		enum Hysteresis {
			HysteresisNothing = 0,		// �q�X�e���V�X�Ȃ��D
			Kakuta = 1,					// �p�c�̎�
		}hysteresis;
		double hysteresis_factor;		// �q�X�e���V�X�����W��

	} tribology;

	struct Cylinder {
		double density;			//�]���̖��x[kg/m^3]
		double young;
		double poisson;
		double ri;
		double ro;
		double l;
		double x0;	// �\���ψʗʁD�i�b�g��p�����ǂ����ɒu�����Ă��������D
		double xg;	// �d�S�̈ʒu�D�i�b�g��p�����ǂ����ɒu�����Ă��������D
		double m;
		double Ix;
		double Iyz;
		struct Spiral {
			double alp;
			double l;
			double r;
			struct Groove {
				double sigma;
				double r;
				double eta[2];	// �aR���S�̂���i0: �Ő����C1:�Đ����j
			}groove[2];
		};
		std::vector<Spiral> spiral;
	} ;
	std::vector<BS_In::Cylinder> nut;
	std::vector<BS_In::Cylinder> shaft;

	struct Circuit {
		struct Ball {
			double density;			//�]���̖��x[kg/m^3]
			double young;
			double poisson;
			double rms;
			double r;
		};
		std::vector<Ball> ball;
		double th0;	// ���׊J�n�ʑ��p
		double th1;	// ���׏I���ʑ��p
		double x0;	// ���׊J�nx���W�iFileIn.cpp ���݂̂Ŏg�p�j
		double x1;	// ���׏I��x���W�iFileIn.cpp ���݂̂Ŏg�p�j
		int is;		// �Ή������ԍ�
		int inut;	// �Ή�����i�b�g�ԍ�
	};

	std::vector<BS_In::Circuit> circuit;

	struct BallCylinderPair {

		struct Groove {
			double mu;
			double zeta;
		}groove[2];
	};
	std::vector<BS_In::BallCylinderPair> BallNutPair;
	std::vector<BS_In::BallCylinderPair> BallShaftPair;

	struct Oil {
		double alpha;
		double beta;
		double k;
		double eta;
		double lm;
	}oil;

	struct Load {
		double x[3];
		double F[3];
		double T[3];
	};
	std::vector<Load> load;

	struct Rigid {
		double g[3];
		double l;
		double t;
	} rigid;

	struct Bound {
		double x0[3];	// [m/s]:	�������͕ψ�
		double ax0[3];	// [rad/s]:	�������͎p���D[1,tan_z,-tan_y]�̍\��
		bool v_const[3];
		bool w_const[3];
	} bound;
};

