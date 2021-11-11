#pragma once
#include <math.h>
#include <Eigen\Dense>
#include "Numeric.h"
#include "Unit.h"
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::VectorXd;
using Eigen::Quaterniond;
using Eigen::Vector4d;


/// <summary>トライボロジー</summary>
/// <remarks>弾性接触や摩擦理論をまとめたクラス．多くは論文の出典が存在する．</remarks>
namespace Tribology {

	void SliceForceRatio(int msmax, double*ForceRatio);
	void SliceForceRatio2d(int nr, int nt, double * r);
	double ReducedYoung(double E0, double m0, double E1, double m1);
	double ReducedMass(double m0, double m1);
	double CompositeRoughness(double s0, double s1);
	double ForceRatio(double lambda);
	double RollerLength(double a, double Rx, double Ry, double xio);

	// クーロン摩擦式．
	class Coulomb {
	public:
		virtual double calc(double mu0, double us, double ur, double s) = 0;
	};
	class CoulombNothing : public Coulomb {
	public:
		double calc(double mu0, double us, double ur, double s);
	};
	class Tangent : public Coulomb {
	public:
		double calc(double mu0, double us, double ur, double s);
	};

	// Hertz の接触計算式．
	class Hertz {
	public:
		virtual void calc(double Rx, double Ry, double dx, double E, double&k, double&a, double&b) = 0;
	};
	class BrewHamrock : public Hertz {
	public:
		void calc(double Rx, double Ry, double dx, double E, double&k, double&a, double&b);
	};

	// 転がり粘性抵抗の式．
	class RollingResistance {
	public:
		virtual double calc(double Rx, double Ry, double bl_D, double a, double b, double w, double E, double um, double fratio, double eta, double alpha, double beta, double k, double xio) = 0;
	};
	class RollingResistanceNothing : public RollingResistance {
	public:
		double calc(double Rx, double Ry, double bl_D, double a, double b, double w, double E, double um, double fratio, double eta, double alpha, double beta, double k, double xio);
	};
	class Houpert : public RollingResistance {
	public:
		double calc(double Rx, double Ry, double bl_D, double a, double b, double w, double E, double um, double fratio, double eta, double alpha, double beta, double k, double xio);
	};
	class AiharaR : public RollingResistance {
	public:
		double calc(double Rx, double Ry, double bl_D, double a, double b, double w, double E, double um, double fratio, double eta, double alpha, double beta, double k, double xio);
	};
	class Fujiwara : public RollingResistance {
	public:
		double calc(double Rx, double Ry, double bl_D, double a, double b, double w, double E, double um, double fratio, double eta, double alpha, double beta, double k, double xio);

	};
	class GoksemAihara : public RollingResistance {
	public:
		double calc(double Rx, double Ry, double bl_D, double a, double b, double w, double E, double um, double fratio, double eta, double alpha, double beta, double k, double xio);
		double Starvation(double xio, double D_GH, double L_GH, double b, double Kt);
	};

	class Hysteresis {
	public:
		virtual double calc(double fh, double b, double Q) = 0;
	};
	class Kakuta : public Hysteresis {
	public:
		double calc(double fh, double b, double Q);
	};
	class HysteresisNothing : public Hysteresis {
	public:
		double calc(double fh, double b, double Q);
	};

	// 油膜厚さの式．
	class FilmThickness {
	public:
		virtual double calc(double w, double E, double Rx, double Ry, double alpha, double eta, double u, double lm, double bh) = 0;
	};
	class FilmThicknessNothing : public FilmThickness {
	public:
		double calc(double w, double E, double Rx, double Ry, double alpha, double eta, double u, double lm, double bh);
	};
	class HamrockDowsonHc : public FilmThickness {
	public:
		double calc(double w, double E, double Rx, double Ry, double alpha, double eta, double u, double lm, double bh);
		double Starvation(double H, double Rx, double lm, double bh);
	};
	class HamrockDowsonHmin : public FilmThickness {
	public:
		double calc(double w, double E, double Rx, double Ry, double alpha, double eta, double u, double lm, double bh);
		double Starvation(double H, double Rx, double lm, double bh);
	};



	// 発熱による補正係数．
	double ErtelGrubin(double eta, double beta, double k, double u);

	// トラクション係数（せん断摩擦係数）の式．
	class Traction {
	public:
		virtual double calc(double eta0, double p0, double v, double u) = 0;
	};
	class AiharaT : public Traction {
	public:
		double calc(double eta0, double p0, double v, double u);
	};

	// 接触に対する減衰力
	class DampingForce {
	public:
		virtual double calc(double k, double zeta, double m, double v, double x) = 0;
	};
	class Tsuji : public DampingForce {
	public:
		double calc(double k, double zeta, double m, double v, double x);
	}; 
	class KelvinVoigt : public DampingForce {
		public:
			double calc(double k, double zeta, double m, double v, double x);
	};


	// スタブ
	class Stab_Traction : public Traction {
	public:
		double calc(double eta0, double p0, double v, double u);
	};
	class Stab_Coulomb : public Coulomb {
	public:
		double calc(double mu0, double us, double ur, double s);
	};
	class Stab_RollingResist : public RollingResistance {
	public:
		double calc(double Rx, double Ry, double bl_D, double a, double b, double w, double E, double um, double fratio, double eta, double alpha, double beta, double k, double xio);
	};
}



//class HamrockDowson_old : public FilmThickness {
//public:
//	double calc(double w, double E, double Rx, double Ry, double alpha, double eta, double u, double lm, double bh);
//};
//class Aihara_old : public RollingResistance {
//public:
//	double calc(double Rx, double Ry, double a, double b, double w, double E, double um, double eta, double alpha, double beta, double k);
//};
