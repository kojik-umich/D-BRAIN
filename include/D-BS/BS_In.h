#pragma once
#include <vector>

class BS_In {
public:
	virtual void Nothing(void) = 0;	// 抽象クラスであることを示すため，意味のない関数を定義する．

public:
	int ballnum;
	struct Tribology {
		enum RollingResistance {
			RollingResistanceNothing = 0,	// 転がり摩擦なし
			GoksemAihara = 1,	// 相原の式
			Aihara = 2,	// 相原の式
			Fijiwara = 3,	// 藤原の式
			Houpert = 4	// Houpertの式
		} rollingresistance;

		enum Coulomb {
			CoulombNothing = 0,		// クーロン摩擦なし．
			Tangent = 1,			// タンジェントカーブ．
			SimpleCoulomb = 2
		} coulomb;
		double coulomb_slope;

		enum Ellipse {
			Mesh1d = 0,	// 1次元メッシュ
			Mesh2d = 1	// 2次元メッシュ
		} ellipse;
		int ellipse_mesh[2];

		enum FilmThickness {
			FilmThicknessNothing = 0,	// 油膜なし．
			HamrockDowsonHc = 1,			// HamrockDowson式による中央油膜厚さ．
			HamrockDowsonHmin = 2,			// HamrockDowson式による最小油膜厚さ．
		} filmThickness;

		enum Hysteresis {
			HysteresisNothing = 0,		// ヒステリシスなし．
			Kakuta = 1,					// 角田の式
		}hysteresis;
		double hysteresis_factor;		// ヒステリシス損失係数

	} tribology;

	struct Cylinder {
		double density;			//転動体密度[kg/m^3]
		double young;
		double poisson;
		double ri;
		double ro;
		double l;
		double x0;	// 予圧変位量．ナット専用だけどここに置かせてください．
		double xg;	// 重心の位置．ナット専用だけどここに置かせてください．
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
				double eta[2];	// 溝R中心のずれ（0: η成分，1:ζ成分）
			}groove[2];
		};
		std::vector<Spiral> spiral;
	} ;
	std::vector<BS_In::Cylinder> nut;
	std::vector<BS_In::Cylinder> shaft;

	struct Circuit {
		struct Ball {
			double density;			//転動体密度[kg/m^3]
			double young;
			double poisson;
			double rms;
			double r;
		};
		std::vector<Ball> ball;
		double th0;	// 負荷開始位相角
		double th1;	// 負荷終了位相角
		double x0;	// 負荷開始x座標（FileIn.cpp 内のみで使用）
		double x1;	// 負荷終了x座標（FileIn.cpp 内のみで使用）
		int is;		// 対応する条番号
		int inut;	// 対応するナット番号
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
		double x0[3];	// [m/s]:	初期入力変位
		double ax0[3];	// [rad/s]:	初期入力姿勢．[1,tan_z,-tan_y]の構成
		bool v_const[3];
		bool w_const[3];
	} bound;
};

