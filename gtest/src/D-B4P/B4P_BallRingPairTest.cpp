#include "pch.h"
#include "B4P_StabIn.h"
#include "B4P_StabOut.h"

class B4P_BallRingPairTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	B4P_BallInnerRingPair BIP;
	B4P_BallOuterRingPair BOP;
	Ball BL;
	B4P_InnerRing IR;
	B4P_OuterRing OR;
	B4P_StabIn FI;
	B4P_StabOut FO;

	// ²σ³π50BSWZ02ΜΰΜΙέθ
	void init_50BSWZ() {
		// 
		FI.LB.eta0 = 0.5;	// SxiPa*sj(K)
		FI.LB.beta0 = 0.5;	// ·xSxWiKj
		FI.LB.k0 = 0.145;	// ϋM`±¦[W/(mEK)]
		FI.LB.alpha0 = 0.2;	// ³ΝSxW[mm2/kgf]
		FI.LB.lm0 = -1;	// jXJX·³


		// Κ¨«l
		double D = 0.0111;				// Κa[m]
		double E = 208000000000;		// O¦[Pa]
		double por = 0.29;				// |A\δ[-]
		double den = 7830;				// §x[kg/m^3]
		double rms = 0.00000002;		// e³rms[m]
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}

		this->BL.init(D, E, por, den, rms, x_const, Rx_const);
		int msmax = 21;
		FO.allocate(1, _MAX_CONTACT_, msmax);

		// ΰΦ¨«l
		FI.ballpcd = 0.0675;
		double clr = 0.000008;
		FI.IR.R[0] = 0.00589;
		FI.IR.R[1] = 0.00589;
		FI.IR.Rox[0] = 0.000167;
		FI.IR.Rox[1] = -0.000167;
		// aRSPCD[m]
		FI.IR.Rod[0] = FI.ballpcd - clr * 0.5 + 2 * sqrt((FI.IR.R[0] - D * 0.5)*(FI.IR.R[0] - D * 0.5) - FI.IR.Rox[0] * FI.IR.Rox[0]);
		FI.IR.Rod[1] = FI.ballpcd - clr * 0.5 + 2 * sqrt((FI.IR.R[1] - D * 0.5)*(FI.IR.R[1] - D * 0.5) - FI.IR.Rox[1] * FI.IR.Rox[1]);

		FI.IR.hedge[0];
		FI.IR.hedge[1];
		FI.IR.rms = 0.00000006;
		FI.IR.m = 1;
		FI.IR.E = 208000000000;
		FI.IR.por = 0.29;
		FI.IR.Ix, FI.IR.Iyz, FI.IR.Iyz;
		FI.BIP.mu = 0.1;
		FI.BIP.dzeta = 0.2;
		// _
		FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::RollingResistanceNothing;
		FI.TB.coulomb = B4P_In::Tribology::Tangent;
		FI.TB.coulomb_slope = 1e100;
		FI.TB.hysteresis = B4P_In::Tribology::HysteresisNothing;
		this->IR.init(FI.IR, FI.ballpcd);
		this->BIP.link(&this->BL, &this->IR);
		this->BIP.init(FI);

		// OΦ¨«l
		FI.OR.Rox[0] = 0.0001665;
		FI.OR.Rox[1] = -0.0001665;
		FI.OR.R[0] = 0.00589;
		FI.OR.R[1] = 0.00589;
		// aRSPCD[m]
		FI.OR.Rod[0] = FI.ballpcd + clr * 0.5 - 2 * sqrt((FI.OR.R[0] - D * 0.5)*(FI.OR.R[0] - D * 0.5) - FI.OR.Rox[0] * FI.OR.Rox[0]);
		FI.OR.Rod[1] = FI.ballpcd + clr * 0.5 - 2 * sqrt((FI.OR.R[1] - D * 0.5)*(FI.OR.R[1] - D * 0.5) - FI.OR.Rox[1] * FI.OR.Rox[1]);
		FI.OR.hedge[0];
		FI.OR.hedge[1];
		FI.OR.rms = 0.00000006;
		FI.OR.m = 1;
		FI.OR.E = 207760000000;
		FI.OR.por = 0.29;
		FI.OR.Ix, FI.OR.Iyz, FI.OR.Iyz;
		FI.LB.eta0 = 0.5;	// SxiPa*sj(K)
		FI.LB.beta0 = 0.5;	// ·xSxWiKj
		FI.LB.k0 = 0.145;	// ϋM`±¦[W/(mEK)]
		FI.LB.alpha0 = 0.2;	// ³ΝSxW[mm2/kgf]
		FI.LB.lm0 = -1;	// jXJX·³
		FI.BOP.mu = 0.1;	// CW
		FI.BOP.dzeta = 0.2;	// ΈW
		this->OR.init(FI.OR, FI.ballpcd);
		this->BOP.link(&this->BL, &this->OR);
		this->BOP.init(FI);


		return;
	};

	// ²σ³π25BSWZ01ΜΰΜΙέθ
	void init_25BSWZ() {
		// όΝf[^ΝΊLΜΰΜπp’½
		// \\ans00978\kiken\BRAIN\\[XR[hΌ\DBRAINpp\docs\06_D4Bv100\05_PΜeXg\B4P_BallRingPair
		// ²σ^Τ@25BSWZ01
		// PairNX¨«l
		this->FI.BOP.dzeta = 0.2;	// Έ
		this->FI.BOP.mu = 0.1;		// CW
		this->FI.BIP.dzeta = 0.2;	// Έ
		this->FI.BIP.mu = 0.1;		// CW

		// Κ¨«l
		double D = 0.00635;				// Κa[m]
		double E = 207760000000;		// O¦[Pa]
		double por = 0.29;				// |A\δ[-]
		double den = 7830;				// §x[kg/m^3]
		double rms = 0.00000002;		// e³rms[m]
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}
		this->BL.init(D, E, por, den, rms, x_const, Rx_const);
		int msmax = 21;
		FO.allocate(1, _MAX_CONTACT_, msmax);

		// ΰΦ¨«l
		FI.ballpcd = 0.0355;
		double clr = 0;
		FI.IR.Rox[0] = 0.0000725;
		FI.IR.Rox[1] = -0.0000725;
		FI.IR.R[0] = 0.00332105;
		FI.IR.R[1] = 0.00332105;
		// aRSPCD[m]iapcdΌaC0.0357535 m ­η’j
		FI.IR.Rod[0] = FI.ballpcd - clr * 0.5 + 2 * sqrt((FI.IR.R[0] - D * 0.5)*(FI.IR.R[0] - D * 0.5) - FI.IR.Rox[0] * FI.IR.Rox[0]);
		FI.IR.Rod[1] = FI.ballpcd - clr * 0.5 + 2 * sqrt((FI.IR.R[1] - D * 0.5)*(FI.IR.R[1] - D * 0.5) - FI.IR.Rox[1] * FI.IR.Rox[1]);
		FI.IR.hedge[0];
		FI.IR.hedge[1];
		FI.IR.rms = 0.00000006;
		FI.IR.m = 1;
		FI.IR.E = 207760000000;
		FI.IR.por = 0.29;
		FI.IR.Ix, FI.IR.Iyz, FI.IR.Iyz;

		// _
		FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::RollingResistanceNothing;
		FI.TB.coulomb = B4P_In::Tribology::Tangent;
		FI.TB.coulomb_slope = 1000;
		FI.TB.hysteresis = B4P_In::Tribology::HysteresisNothing;
		this->IR.init(FI.IR, FI.ballpcd);
		this->BIP.link(&this->BL, &this->IR);
		this->BIP.init(FI);

		// OΦ¨«l
		FI.OR.Rox[0] = 0.0000725;
		FI.OR.Rox[1] = -0.0000725;
		FI.OR.R[0] = 0.00332105;
		FI.OR.R[1] = 0.00332105;
		FI.OR.Rod[0] = FI.ballpcd + clr * 0.5 - 2 * sqrt((FI.OR.R[0] - D * 0.5)*(FI.OR.R[0] - D * 0.5) - FI.OR.Rox[0] * FI.OR.Rox[0]);
		FI.OR.Rod[1] = FI.ballpcd + clr * 0.5 - 2 * sqrt((FI.OR.R[1] - D * 0.5)*(FI.OR.R[1] - D * 0.5) - FI.OR.Rox[1] * FI.OR.Rox[1]);
		FI.OR.hedge[0];
		FI.OR.hedge[1];
		FI.OR.rms = 0.00000006;
		FI.OR.m = 1;
		FI.OR.E = 207760000000;
		FI.OR.por = 0.29;
		FI.OR.Ix, FI.OR.Iyz, FI.OR.Iyz;
		FI.LB.eta0 = 0.5;	// SxiPa*sj(K)
		FI.LB.beta0 = 0.5;	// ·xSxWiKj
		FI.LB.k0 = 0.145;	// ϋM`±¦[W/(mEK)]
		FI.LB.alpha0 = 0.2;	// ³ΝSxW[mm2/kgf]
		FI.LB.lm0 = -1;	// jXJX·³
		this->OR.init(FI.OR, FI.ballpcd);
		this->BOP.link(&this->BL, &this->OR);
		this->BOP.init(FI);

	};

	virtual void SetUp() {
		Vector3d x = Vector3d(0, 0, 0);
		Vector3d v = Vector3d(0, 0, 0);
		Quaterniond q = Quaterniond(1, 0, 0, 0); // ϊNH[^jIΝSΔ [0,0,0,1] Ζ·ιDi½Ύ΅RXgN^ΜdlγΤͺΩΘιj
		Vector3d w = Vector3d(0, 0, 0);
		FI.msmax = 21;
		Rigid::l = 1;
		Rigid::t = 1;
		Rigid::g = Vector3d(0, 0, 0);
		this->OR.set_y(x, v, q, w);
		this->IR.set_y(x, v, q, w);
		int msmax = 21;
		FO.allocate(1, _MAX_CONTACT_, msmax);
	}

	virtual void TearDown() {

	}
};



// (1) ΫνEΚEΰOΦ©]¬xπp[^Ζ΅½?πμ¬΅CδrD
TEST_F(B4P_BallRingPairTest, get_us_ur_1) {
	this->init_25BSWZ();
	double wc = 44.4604328912874;		// Ϋνφ]¬x[rad/s]
	double wrn_0 = 0.0;					// OΦρ]¬x[rad/s]
	double wrn_1 = 104.719755119660;	// ΰΦρ]¬x[rad/s]
	double omega_x = -223.906546549487;	// Κ©]¬x Φx[rad/s]
	double omega_z = 287.884790174018;  // Κ©]¬x Φz[rad/s]
	double dx_out = 4.389357066845456E-006;		// OΦ€ΪίΚ(balacvZΚ)
	double dx_in = 4.518845377961239E-006;		// ΰΦ€ΪίΚ(balacvZΚ)
	double sin_alp1 = 0.537033675917149;		// OΦΪGpΜ³·(balacvZΚ)
	double sin_alp2 = 0.537125019737040;		// ΰΦΪGpΜ³·(balacvZΚ)
	double cos_alp1 = 0.843560804525029;		// OΦΪGpΜ]·(balacvZΚ)
	double cos_alp2 = -0.843502645622694;		// ΰΦΪGpΜ]·(balacvZΚ)

	// ΚΞΚuΝΪίΚ©ηtZ
	double R_out = this->OR.GV[0].r - this->BL.r;
	double R_in = this->IR.GV[0].r - this->BL.r;
	Vector3d er = Vector3d(0, 0, 1);			// aϋόxNg
	Vector3d bl_x_in = Vector3d(0, 0, 0);
	Vector3d bl_x_out = Vector3d(0, 0, 0);

	Vector3d x = Vector3d(0, 0, 0);
	Vector3d v = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // ϊNH[^jIΝSΔ [0,0,0,1] Ζ·ιDi½Ύ΅RXgN^ΜdlγΤͺΩΘιj
	Vector3d wo = Vector3d(wrn_0, 0, 0);
	this->OR.set_y(x, v, q, wo);
	Vector3d wi = Vector3d(wrn_1, 0, 0);
	this->IR.set_y(x, v, q, wi);
	Vector3d bl_x = Vector3d(0, 0, 0.5 * FI.ballpcd);
	Vector3d bl_v = Vector3d(0, -0.5 * FI.ballpcd * wc, 0);
	Vector3d bl_w = Vector3d(omega_x, 0, omega_z);
	this->BL.set_y(bl_x, bl_v, q, bl_w);

	// ΚSΙΞ·ιΪG_ΜΞ£
	Vector3d _p_out = Vector3d(this->BL.r * sin_alp1, 0, this->BL.r * cos_alp1 + 0.5 * FI.ballpcd);
	Vector3d _p_in = Vector3d(this->BL.r * sin_alp2, 0, this->BL.r * cos_alp2 + 0.5 * FI.ballpcd);
	Vector3d uso_, uro_, usi_, uri_;
	this->BOP.get_us_ur(_p_out, uso_, uro_);

	double uso_norm = uso_.norm(), uro_norm = uro_.norm();
	this->BIP.get_us_ur(_p_in, usi_, uri_);
	double usi_norm = usi_.norm(), uri_norm = uri_.norm();

	// baltacΙΐ³κΔ’½?πQlΙ©μ΅½?(κΟΜθ`ϋ@ͺΩΘι½ίόΟ)
	// iΪG_Μ¬xΜvZϋ@ΰ{vOΜvZϋ@ΙνΉ±ρΕ’ιj
	double u1, u2, baltac_uso, baltac_usi, baltac_Uo = 0, baltac_Ui = 0;
	u1 = 0.5 * FI.ballpcd  * wc - 0.5 * (FI.ballpcd + this->BL.D * cos_alp1) * wrn_0;
	u2 = 0.5 * this->BL.D *(-omega_x * cos_alp1 + omega_z * sin_alp1);
	baltac_uso = u1 - u2;
	baltac_Uo = (u1 + u2) * 0.5;
	u1 = 0.5 * FI.ballpcd  * wc - 0.5 * (FI.ballpcd + this->BL.D * cos_alp2) * wrn_1;
	u2 = 0.5 * this->BL.D *(-omega_x * cos_alp2 + omega_z * sin_alp2);
	baltac_usi = u1 - u2;
	baltac_Ui = (u1 + u2) * 0.5;

	double err = 0.05;

	cout << uro_norm << "\t" << abs(baltac_Uo) << endl;
	cout << uso_norm << "\t" << abs(baltac_uso) << endl;
	cout << uri_norm << "\t" << abs(baltac_Ui) << endl;
	cout << usi_norm << "\t" << abs(baltac_usi) << endl;
	// eλ·
	EXPECT_NEAR(uro_norm, abs(baltac_Uo), abs(baltac_Uo)*err);
	EXPECT_NEAR(uso_norm, abs(baltac_uso), abs(baltac_uso)*err);
	EXPECT_NEAR(uri_norm, abs(baltac_Ui), abs(baltac_Ui)*err);
	EXPECT_NEAR(usi_norm, abs(baltac_usi), abs(baltac_usi)*err);

	return;
}

// (1) ΫνEΚEΰOΦ©]¬xπp[^Ζ΅½?πμ¬΅CδrD
TEST_F(B4P_BallRingPairTest, get_us_ur2_1) {

	this->init_25BSWZ();
	double wc = 44.4604328912874;		// Ϋνφ]¬x[rad/s]
	double wrn_0 = 0.0;					// OΦρ]¬x[rad/s]
	double wrn_1 = 104.719755119660;	// ΰΦρ]¬x[rad/s]

	double omega_x = -223.906546549487;	// Κ©]¬x Φx[rad/s]
	double omega_z = 287.884790174018;  // Κ©]¬x Φz[rad/s]
	double dx_out = 4.389357066845456E-006;		// OΦ€ΪίΚ(balacvZΚ)
	double dx_in = 4.518845377961239E-006;		// ΰΦ€ΪίΚ(balacvZΚ)
	double sin_alp1 = 0.537033675917149;		// OΦΪGpΜ³·(balacvZΚ)
	double sin_alp2 = 0.537125019737040;		// ΰΦΪGpΜ³·(balacvZΚ)
	double cos_alp1 = 0.843560804525029;		// OΦΪGpΜ]·(balacvZΚ)
	double cos_alp2 = -0.843502645622694;		// ΰΦΪGpΜ]·(balacvZΚ)

	// ΚΞΚuΝΪίΚ©ηtZ
	double R_out = this->OR.GV[0].r - this->BL.r;
	double R_in = this->IR.GV[0].r - this->BL.r;
	Vector3d er = Vector3d(0, 0, 1);			// aϋόxNg
	Vector3d bl_x_in = Vector3d(0, 0, 0.5 * FI.ballpcd);
	Vector3d bl_x_out = Vector3d(0, 0, 0.5 * FI.ballpcd);

	Vector3d x = Vector3d(0, 0, 0);
	Vector3d v = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // ϊNH[^jIΝSΔ [0,0,0,1] Ζ·ιDi½Ύ΅RXgN^ΜdlγΤͺΩΘιj
	Vector3d wo = Vector3d(wrn_0, 0, 0);
	this->OR.set_y(x, v, q, wo);
	Vector3d wi = Vector3d(wrn_1, 0, 0);
	this->IR.set_y(x, v, q, wi);
	Vector3d bl_v = Vector3d(0, -0.5 * FI.ballpcd * wc, 0);
	Vector3d bl_w = Vector3d(omega_x, 0, omega_z);
	this->BL.set_y(bl_x_in, bl_v, q, bl_w);

	// ΚSΙΞ·ιΪG_ΜΞ£
	Vector3d _p_out = Vector3d(this->BL.r * sin_alp1, 0, this->BL.r * cos_alp1 + 0.5 * FI.ballpcd);
	Vector3d _p_in = Vector3d(this->BL.r * sin_alp2, 0, this->BL.r * cos_alp2 + 0.5 * FI.ballpcd);
	Vector3d uso_, uro_, usi_, uri_;
	this->BOP.get_us_ur2(_p_out, er, uso_, uro_);

	double uso_norm = uso_.norm(), uro_norm = uro_.norm();
	this->BIP.get_us_ur2(_p_in, er, usi_, uri_);
	double usi_norm = usi_.norm(), uri_norm = uri_.norm();

	// baltacΙΐ³κΔ’½?πQlΙ©μ΅½?(κΟΜθ`ϋ@ͺΩΘι½ίόΟ)
	// iΪG_Μ¬xΜvZϋ@ΰ{vOΜvZϋ@ΙνΉ±ρΕ’ιj
	double u1, u2, baltac_uso, baltac_usi, baltac_Uo = 0, baltac_Ui = 0;
	u1 = 0.5 * (FI.ballpcd + this->BL.D * cos_alp1)*(wc - wrn_0);
	u2 = 0.5 * this->BL.D *(-(omega_x - (wc - wrn_0))* cos_alp1 + omega_z * sin_alp1);
	baltac_uso = u1 - u2;
	baltac_Uo = (u1 + u2) * 0.5;
	u1 = 0.5 * (FI.ballpcd + this->BL.D * cos_alp2)*(wc - wrn_1);
	u2 = 0.5 * this->BL.D *(-(omega_x - (wc))* cos_alp2 + omega_z * sin_alp2);
	baltac_usi = u1 - u2;
	baltac_Ui = (u1 + u2) * 0.5;

	double err = 0.05;

	cout << uro_norm << "\t" << abs(baltac_Uo) << endl;
	cout << uso_norm << "\t" << abs(baltac_uso) << endl;
	cout << uri_norm << "\t" << abs(baltac_Ui) << endl;
	cout << usi_norm << "\t" << abs(baltac_usi) << endl;
	// eλ·
	EXPECT_NEAR(uro_norm, abs(baltac_Uo), abs(baltac_Uo)*err);
	EXPECT_NEAR(uso_norm, abs(baltac_uso), abs(baltac_uso)*err);
	EXPECT_NEAR(uri_norm, abs(baltac_Ui), abs(baltac_Ui)*err);
	EXPECT_NEAR(usi_norm, abs(baltac_usi), abs(baltac_usi)*err);

	return;
}


// (1) CWπΕθlΙ΅½Ζ«ΙΚ[ΰOΦΤΜΧdͺθvZΖ―ΆΙΘΑΔ’ι©Ψ
TEST_F(B4P_BallRingPairTest, calc_force_1) {
	this->init_25BSWZ();

	// ΪGσΤ
	double sin_alp1 = 0.540517404112816;		// OΦΪGpΜ³·(balacvZΚ)
	double cos_alp1 = 0.841332832980588;		// OΦΪGpΜ]·(balacvZΚ)
	double dx_out = 4.389357066845456E-006;		// OΦ€ΪίΚ(balacvZΚ)
	double _Qo = 21.3811947318042;				// ]?ΜΧd[kgf](balacvZΚ)

	this->BOP.TR = new Tribology::Stab_Traction();	// gNVWπ0.2ΙΕθ
	this->BOP.CL = new Tribology::Stab_Coulomb();	// N[CWπ0.3ΙΕθ
	this->BOP.RR = new Tribology::Stab_RollingResist();	// ]ͺθS«οRπθΙΕθ
	this->BOP.FT = new Tribology::FilmThicknessNothing();	// ϋϊ³π0ΙΕθ

	// Κͺ³ρ]Εΐi^?΅Δ’ικπΌθ
	int i = 1;			// ΪG·ιaΤ
	Vector3d x = Vector3d(0, 0, 0), v = Vector3d(0, 0, 0), w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	this->OR.set_y(x, v, q, w);
	Vector3d er = Vector3d(0, 0, 1);
	double dx = this->OR.GV[i].Rx;
	double R_out = this->OR.GV[i].r + dx_out - this->BL.r;			// aRS©ηΚSΜ£
	double rp = this->OR.GV[i].r + dx_out * 0.5;					// aRS©ηΪG_Μ£
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + dx, 0, cos_alp1 * R_out) + er * this->OR.GV[i].Rr;
	Vector3d bl_v = Vector3d(0, -1, 0);
	this->BL.set_y(bl_x_out, bl_v, q, w);


	Vector3d Fbi, Tbi, Fib, Tib;

	this->BOP.calc_force(Fbi, Tbi, Fib, Tib);

	// θvZΜΚΖδr
	double _fn = _Qo * 9.8;
	Vector3d _Fn = Vector3d(-sin_alp1, 0, -cos_alp1) * _fn;	// ΌΧd[N]
	Vector3d _Fs = Vector3d(0, _fn * 0.3, 0);				// θC[N]
	Vector3d _Fbi = _Fn + _Fs;
	double err = 0.05;
	EXPECT_NEAR((Fbi - _Fbi).norm(), 0, _Fbi.norm() * err);
	double _rb = this->BL.r - dx_out * 0.5;
	double _Ts_norm = _fn * 0.3 * _rb;
	Vector3d _Ts = Vector3d(-cos_alp1, 0, sin_alp1) * _Ts_norm;			// θC[N]
	Vector3d _Tr = Vector3d(cos_alp1, 0, -sin_alp1) * 0.1;
	Vector3d _Tbi = _Ts + _Tr;
	EXPECT_NEAR((Tbi - _Tbi).norm(), 0, _Tbi.norm() * err);

	// write()ΜeXgπΛΔΪGΘ~βΪGpπΨ

	this->BOP.write(FO.BOP[0]);
	Vector3d Fn(FO.BOP[0].GV[i].Fn);		// ΌΧd[N]
	EXPECT_NEAR((Fn - _Fn).norm(), 0, _Fn.norm() * err);
	Vector3d Fs(FO.BOP[0].GV[i].Fs);		// θC[N]
	EXPECT_NEAR((Fs - _Fs).norm(), 0, _Fs.norm() * err);
	double ea = FO.BOP[0].GV[i].a;											// Θ~·Όa[mm](ΐΫ)
	double _ea = 0.668610936573505e-3;										// Θ~·Όa[mm](ϊ?l)
	EXPECT_NEAR(_ea, ea, _ea * err);
	double eb = FO.BOP[0].GV[i].b;											// Θ~ZΌa[mm](ΐΫ)
	double _eb = 9.720461846754701e-5;								// Θ~ZΌa[mm](ϊ?l)
	EXPECT_NEAR(_eb, eb, _eb * err);
	double alp = FO.BOP[0].GV[i].phi;											// ΪGp[](ΐΫ)
	double _alp = Unit::deg2rad(32.7188677338554);										// ΪGp[](ϊ?l)
	EXPECT_NEAR(_alp, alp, _alp * err);

	return;
}

// (1) ΪGp0Ι¨―ιΚΪίΚΖeνϋόxNgπmF
TEST_F(B4P_BallRingPairTest, how_Contact_1) {
	this->init_25BSWZ();

	// ΪGp0EΪίΚ1mmΜΖ«ΜΚSΚuπtZ
	// ΚSyΐW = ΪίΚ + aRΌa - ΚΌa + aRSyΐW
	// ½Ύ΅CaSΚuΝ(0.0000725, 0.035246430660370774* 0.5, 0)
	// ΪG·ιaΝ +X €Ζ·ι

	Vector3d x(-0.0000725, (0.00332105 + 0.035246430660370774 * 0.5 - 0.00635 * 0.5 + 0.001), 0);

	Vector3d er, eg;
	double dx;

	int i = 1; 
	bool c = this->BOP.how_Contact(i, x, er, eg, dx);

	double err = 0.01;					// eλ·
	Vector3d expcted_er(0, 1, 0);		// ΚΚxNg
	Vector3d expcted_eg(0, 1, 0);		// ΚΪGϋόxNg
	double expected_dx = 0.001;			// ΚΪίΚ

	EXPECT_NEAR((er - expcted_er).norm(), 0, er.norm()*err);
	EXPECT_NEAR((eg - expcted_eg).norm(), 0, eg.norm()*err);
	EXPECT_NEAR(dx, expected_dx, dx*err);

	return;
};

// (2) {[Sͺg[X©ηoΔ’½ηΪGΘ΅Μ»θΖΘι±ΖπmF
TEST_F(B4P_BallRingPairTest, how_Contact_2) {
	this->init_25BSWZ();

	// {[Sπg[XζθO€ΙέθCΪG·ιaΝ+X€Ιέθ
	Vector3d x(-0.0000725, 1.0, 0);
	int i = 1;

	Vector3d er, eg;
	double dx;
	bool c = this->BOP.how_Contact(i, x, er, eg, dx);

	EXPECT_EQ(false, c);
	return;
};

// (3) H’έΚͺlΜκΪGΘ΅Μ»θΖΘι±ΖπmF
TEST_F(B4P_BallRingPairTest, how_Contact_3) {
	this->init_25BSWZ();

	// ΪίΚ-10um©ΒΪGp0ΖΘιζ€ΙΚΚuπwθ
	Vector3d x(-0.0000725, (0.00332105 + 0.035246430660370774 * 0.5 - 0.00635 * 0.5 - 0.0001), 0);
	//Vector3d v(0.0, -0.01, 0.0);
	int i = 1;

	Vector3d er, eg;
	double dx;

	bool c = this->BOP.how_Contact(i, x,  er, eg, dx);
	EXPECT_EQ(false, c);
	return;
};


// (4) ΚͺaS_ΖdΘιΖ«CΪGΘ΅ΙΘι©mF
TEST_F(B4P_BallRingPairTest, how_Contact_4) {
	this->init_25BSWZ();

	// ΚΚuπaSΚuΙzu
	Vector3d x(-0.0000725, 0.035246430660370774 * 0.5, 0);
	int i = 1;

	Vector3d er, eg;
	double dx;

	// ΚνΝ θ¦Θ’ͺCΚaπ1.2{Ι΅ΔΚͺaSΙ ιΕΰOΦΖΪG·ιζ€Ι²?
	double D = 0.00635 * 1.2;		// Κa[m]
	double E = 207760000000;		// O¦[Pa]
	double por = 0.29;				// |A\δ[-]
	double den = 7830;				// §x[kg/m^3]
	double rms = 0.00000002;		// e³rms[m]
	bool x_const[3], Rx_const[3];
	for (int i = 0; i < 3; i++) {
		x_const[i] = false;
		Rx_const[i] = false;
	}
	this->BL.init(D, E, por, den, rms, x_const, Rx_const);


	bool c = this->BOP.how_Contact(i, x, er, eg, dx);

	double err = 0.05; // eλ·
	Vector3d expcted_er(0, 1, 0);
	Vector3d expcted_eg(0, 1, 0);
	double expected_dx = 0.001;

	EXPECT_EQ(false, c);
	return;
};

// (5) ΚpΖΪGpͺρ[©ΒCΪGEρΪGΜ«EtίΜ?πmF
TEST_F(B4P_BallRingPairTest, how_Contact_5) {
	this->init_25BSWZ();

	// ΚπΚp60CΪGp30CΪίΚ0ΜΚuΙzu
	// xΐW = aRSxΐW + (ΪίΚ + aRΌa - ΚΌa) * sin30
	// yΐW = ((ΪίΚ + aRΌa - ΚΌa) * cos30 + aPCDΌa) * sin60
	// zΐW = ((ΪίΚ + aRΌa - ΚΌa) * cos30 + aPCDΌa) * cos60

	Vector3d x(5.25E-07, 0.01537169, 0.008874849);
	int i = 1;
	Vector3d er, eg;
	double dx;

	// ΪίΚ0ΜΚu©η +x, +y, +z ϋόΙχ¬ΚΪ?³ΉιΖΪG·ι©»θ
	Vector3d ep_x(1e-6, 0, 0), ep_y(0, 1e-6, 0), ep_z(0, 0, 1e-6);
	bool cx1 = this->BOP.how_Contact(i, x + ep_x, er, eg, dx);
	bool cy1 = this->BOP.how_Contact(i, x + ep_y, er, eg, dx);
	bool cz1 = this->BOP.how_Contact(i, x + ep_z, er, eg, dx);
	EXPECT_EQ(true, cx1);
	EXPECT_EQ(true, cy1);
	EXPECT_EQ(true, cz1);

	// ΪίΚ0ΜΚu©η -x, -y, -z ϋόΙχ¬ΚΪ?³ΉιΖρΪGΙΘι©»θ
	bool cx2 = this->BOP.how_Contact(i, x - ep_x, er, eg, dx);
	bool cy2 = this->BOP.how_Contact(i, x - ep_y, er, eg, dx);
	bool cz2 = this->BOP.how_Contact(i, x - ep_z, er, eg, dx);
	EXPECT_EQ(false, cx2);
	EXPECT_EQ(false, cy2);
	EXPECT_EQ(false, cz2);
	return;
};


// (1) Κ¬xπΟ¦ΔvZπs’CΈπlΆ΅Θ’κΖΜδrπs€D
TEST_F(B4P_BallRingPairTest, calc_DynamicHertz_1) {
	this->init_25BSWZ();


	/*****d4bόΝπ(²σ³ΝSetUpΦΕέθ)*****/
	// Κͺ²ϋό©η©Δ +y ²ϋόΙ ικπzθD
	// Κ¬xͺ +y ϋόiOΦO€Ιisj

	int i = 0;									// aΤ
	double sin_alp1 = 0.540517404112816;		// OΦΪGpΜ³·(balacvZΚ)
	double cos_alp1 = 0.841332832980588;		// OΦΪGpΜ]·(balacvZΚ)
	Vector3d er= Vector3d(0, 1, 0);					// aϋόxNg
	Vector3d eg_out= Vector3d(sin_alp1, cos_alp1, 0);	// aS¨{[ϋόxNg(baltacΪGp©ηtZ)
	double dx_out= 4.389357066845456E-006;		// OΦ€ΪίΚ(balacvZΚ)
	double Rx_out, Ry_out, cos_alp, sin_alp, a_out, b_out;
	Vector3d p_out;
	// ΚΞΚuΝΪίΚ©ηtZ
	// ΚπΚp 90 (+yϋό)CΪGp 30CΪίΚ 4.389 um ΜΚuΙzu
	// xΐW = aRSxΐW + (ΪίΚ + aRΌa - ΚΌa) * sin30
	// yΐW = ((ΪίΚ + aRΌa - ΚΌa) * cos30 + aPCDΌa) * sin90
	// zΐW = ((ΪίΚ + aRΌa - ΚΌa) * cos30 + aPCDΌa) * cos90
	double R_out = this->OR.GV[0].r + dx_out - this->BL.r;
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + this->OR.GV[0].Rx, cos_alp1 * R_out, 0) + er * this->OR.GV[0].Rr;
	
	// δrΞΫΖΘιΈΘ΅ΜκπζΙvZD
	double k_out, k_in;
	this->BOP.how_Contact(i, bl_x_out, er, eg_out, dx_out);
	this->BOP.calc_Hertz(i, bl_x_out, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out, k_out);
	double _Qout = k_out * pow(dx_out, 1.5);
	// (1) ΚͺΪGΚΙίΓ­ό«Ιis·ικCΪGΝͺε«­Θι©mF
	Vector3d bl_v1 = Vector3d(1, 0, 0);		// Κ¬xiΗΚ©ηίΓ­ό«j
	double Qout1 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v1, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_GE(Qout1, _Qout);	

	// (2) ΚͺΪGΚ©η£κιό«Ιis·ικCΪGΝͺ¬³­Θι©mF
	Vector3d bl_v2 = Vector3d(-1, 0, 0);	// Κ¬xiΗΚ©η£κιό«j
	double Qout2 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v2, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_LE(Qout2, _Qout);

	// (3) ΚͺΪGΚΖΐsΘό«Ιis·ικCΪGΝͺΟνηΘ’±ΖπmF(z¬ͺΜέ)
	Vector3d bl_v3 = Vector3d(0, 0, -1);	// Κ¬xiΗΚΐsϋόj
	double Qout3 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v3, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_NEAR(Qout3, _Qout, 1e-6);

	// (4) ΚͺΪGΚΖΐsΘό«Ιis·ικCΪGΝͺΟνηΘ’±ΖπmF(xy¬ͺΜέ)
	Vector3d bl_v4 = Vector3d(0.841332832980588, -0.540517404112816, 0);	// Κ¬xiΗΚΐsϋόj
	double Qout4 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v4, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_NEAR(Qout4, _Qout, 1e-6);

	// (5) ΈΝΖHertzΪGΝΜvlͺΜC0Ιβ³΅Δ’ι±ΖπmF
	Vector3d bl_v5 = Vector3d(-1e9, 0, 0);		// Κ¬xiΗΚ©η£κιό«j
	double Qout5 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v5, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_NEAR(Qout5, 0, 1e-6);	

	return;
};




// (1) Έ¬ͺͺΘ’Ζ«ΜwcΪGΝπvZ΅CbaltacΜvZΚΖδr
// baltacΖ―ΜΪίΚEΪGpπόΝ΅½πΕC]?ΜΧdEΘ~Όaπ]Ώ
// iόΝπΝALVAΧdΜκπp’ιDΈΝlΆ΅Θ’Dj
TEST_F(B4P_BallRingPairTest, calc_Hertz_1) {
	this->init_25BSWZ();

	// ²σ³ΝSetUpΦΕέθ
	// d4bόΝπ
	int i = 0;									// aΤ
	double sin_alp1 = 0.540517404112816;		// OΦΪGpΜ³·(balacvZΚ)
	double sin_alp2 = 0.540517404112816;		// ΰΦΪGpΜ³·(balacvZΚ)
	double cos_alp1 = 0.841332832980588;		// OΦΪGpΜ]·(balacvZΚ)
	double cos_alp2 = -0.841332832980588;		// ΰΦΪGpΜ]·(balacvZΚ)
	Vector3d er = Vector3d(0, 1, 0);			// aϋόxNg
	Vector3d eg_out = Vector3d(sin_alp1, cos_alp1, 0); // aS¨{[ϋόxNg(baltacΪGp©ηtZ)
	Vector3d eg_in = Vector3d(sin_alp2, cos_alp2, 0); // aS¨{[ϋόxNg(baltacΪGp©ηtZ)
	double dx_out = 4.389357066845456E-006;		// OΦ€ΪίΚ(balacvZΚ)
	double dx_in = 4.518845377961239E-006;		// ΰΦ€ΪίΚ(balacvZΚ)
	double Rx_out, Ry_out, Rx_in, Ry_in, cos_alp, sin_alp, a_out, b_out, k_out, a_in, b_in, k_in;
	Vector3d p_out, p_in;
	// ΚΞΚuΝΪίΚ©ηtZ
	double R_out = this->OR.GV[0].r + dx_out - this->BL.r;
	double R_in = this->IR.GV[0].r + dx_in - this->BL.r;
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + this->OR.GV[0].Rx, cos_alp1 * R_out, 0) + er * this->OR.GV[0].Rr;
	Vector3d bl_x_in = Vector3d(sin_alp2 * R_in + this->IR.GV[0].Rx, cos_alp2 * R_in, 0) + er * this->IR.GV[0].Rr;
	// eXgΞΫΦΜΔΡo΅
	this->BOP.calc_Hertz(i, bl_x_out, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out, k_out);
	double Qout = k_out * pow(dx_out, 1.5);
	this->BIP.calc_Hertz(i, bl_x_in, er, eg_in, dx_in, Rx_in, Ry_in, p_in, cos_alp, sin_alp, a_in, b_in, k_in);
	double Qin = k_in * pow(dx_in, 1.5);

	// δrΞΫΖΘιbalacvZΚifobO[hΙΔζΎj
	// 1: OΦC2:ΰΦ
	double ea1 = 0.668610936573505;				// Θ~·Όa[mm]
	double ea2 = 0.685074685553727;				// Θ~·Όa[mm]
	double eb1 = 9.720461846754701E-002;		// Θ~ZΌa[mm]
	double eb2 = 8.272063984522840E-002;		// Θ~ZΌa[mm]
	double _Q1 = 21.3811947318042;				// ]?ΜΧd[kgf]
	double _Q2 = 21.3811947318041;				// ]?ΜΧd[kgf]
	double _rho_0 = 2 / 6.35;					// ]?ΜΘ¦[1/mm]
	double _rho3_in = -0.301109588834856;		// aΘ¦[1/mm]
	double _rho2_in = 5.579585936574204E-002;	// OΉΜΘ¦(φ]OΉΪGΌaΜtj[1/mm]
	double _Rx_in = 1 / (_rho_0 + _rho2_in);
	double _Ry_in = 1 / (_rho_0 + _rho3_in);

	double _rho3_out = -0.301109588834856;
	double _rho2_out = -4.119892685701449E-002;
	double _Rx_out = 1 / (_rho_0 + _rho2_out);
	double _Ry_out = 1 / (_rho_0 + _rho3_out);
	// ΪGΚu(eXgP[X)ΝΪίΚΖΪGp©ηvZ
	Vector3d _p_out = Vector3d(sin_alp1 * (this->OR.GV[0].r + dx_out * 0.5) + this->OR.GV[0].Rx, cos_alp1 * (this->OR.GV[0].r + dx_out * 0.5), 0)
		+ er * this->OR.GV[0].Rr;
	Vector3d _p_in = Vector3d(sin_alp2 * (this->IR.GV[0].r + dx_in * 0.5) + this->IR.GV[0].Rx, cos_alp2 * (this->IR.GV[0].r + dx_in * 0.5), 0)
		+ er * this->IR.GV[0].Rr;


	cout << "OΦ€Κ" << endl;
	cout << "]?ΜΧdF" << Qout << "\tbaltacF" << _Q1 * 9.8 << endl;
	cout << "ΪGΘ~·ΌaF" << a_out << "\tbaltacF" << ea1 / 1000 << endl;
	cout << "ΪGΘ~ZΌaF" << b_out << "\tbaltacF" << eb1 / 1000 << endl;
	cout << "Θ¦RxF" << Rx_out << "\tbaltacF" << _Rx_out / 1000 << endl;
	cout << "Θ¦RyF" << Ry_out << "\tbaltacF" << _Ry_out / 1000 << endl;
	cout << "ΪG_Κuλ·F" << (p_out - _p_out).norm() << endl;
	//cout << "p = " << p_out << "\thandF" << _p_out <<endl;
	cout << "ΰΦ€Κ" << endl;
	cout << "]?ΜΧdF" << Qin << "\tbaltacF" << _Q2 * 9.8 << endl;
	cout << "ΪGΘ~·ΌaF" << a_in << "\tbaltacF" << ea2 / 1000 << endl;
	cout << "ΪGΘ~ZΌaF" << b_in << "\tbaltacF" << eb2 / 1000 << endl;
	cout << "Θ¦RxF" << Rx_in << "\tbaltacF" << _Rx_in / 1000 << endl;
	cout << "Θ¦RyF" << Ry_in << "\tbaltacF" << _Ry_in / 1000 << endl;
	cout << "ΪG_Κuλ·F" << (p_in - _p_in).norm() << endl;
	//cout << "p = " << p_in << "\thandF" << _p_in << endl;
	// baltacΕΝBrew-HamrockΜ?ΕΝΘ­C
	// ³ΐWJΙζιίΕίΔ’ι½ί?SΙκv΅Θ’D
	// ibaltacΰΜK(k')ΖE(k')ΜvZΝ^?ΝwUp12ΙfΪj
	// λ·ͺ2%ΘΰΕ κΞΒΖ΅½D
	double err = 0.02;
	EXPECT_NEAR(Qout, _Q1*9.8, Qout*err);
	EXPECT_NEAR(Qin, _Q2*9.8, Qin*err);
	EXPECT_NEAR(a_out, ea1 / 1000, a_out*err);
	EXPECT_NEAR(a_in, ea2 / 1000, a_in*err);
	EXPECT_NEAR(b_out, eb1 / 1000, b_out*err);
	EXPECT_NEAR(b_in, eb2 / 1000, b_in*err);
	EXPECT_NEAR(Rx_out, _Rx_out / 1000, Rx_out*err);
	EXPECT_NEAR(Rx_in, _Rx_in / 1000, Rx_in*err);
	EXPECT_NEAR(Ry_out, _Ry_out / 1000, Ry_out*err);
	EXPECT_NEAR(Ry_in, _Ry_in / 1000, Ry_in*err);
	EXPECT_NEAR((p_out - _p_out).norm(), 0, p_out.norm()*err);
	EXPECT_NEAR((p_in - _p_in).norm(), 0, p_in.norm()*err);
	return;
};

// (1)ΚͺΰOΦ~όϋόΙΪ?΅½Ζ«ΜθCπvZ΅CθvZΜlΖδr
// ½Ύ΅CgNVWEN[CWΝθlΙΕθ
TEST_F(B4P_BallRingPairTest, calc_Sliding_1) {
	this->init_25BSWZ();

	int i = 0;	// aΤ
	double a = 0.668610936573505e-3; // ΪGΘ~·a[m]
	double Pmax = 1;		// ΕεΚ³igνΘ’j
	double F_norm = 100;	// ]?ΜΧd[N]
	double fratio = 0.4;	// ϋΪG
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double sin_alp1 = 0.540517404112816;		// OΦΪGpΜ³·(balacvZΚ)
	double sin_alp2 = 0.540517404112816;		// ΰΦΪGpΜ³·(balacvZΚ)
	double cos_alp1 = 0.841332832980588;		// OΦΪGpΜ]·(balacvZΚ)
	double cos_alp2 = -0.841332832980588;		// ΰΦΪGpΜ]·(balacvZΚ)
	double dx_out = 4.389357066845456E-006;		// OΦ€ΪίΚ(balacvZΚ)


	this->BOP.TR = new Tribology::Stab_Traction();	// gNVWπ0.2ΙΕθ
	this->BOP.CL = new Tribology::Stab_Coulomb();	// N[CWπ0.3ΙΕθ

	// ΚͺYϋόiόϋόjΙ³ρ]ΕΑΔ’ισΤπΌθ
	Vector3d x = Vector3d(0, 0, 0);
	Vector3d v = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d w = Vector3d(0, 0, 0);
	this->OR.set_y(x, v, q, w);
	Vector3d bl_v = Vector3d(0, 1, 0);
	Vector3d er = Vector3d(0, 0, 1);
	double R_out = this->OR.GV[0].r + dx_out - this->BL.r;
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + this->OR.GV[0].Rx, 0, cos_alp1 * R_out) + er * this->OR.GV[0].Rr;
	Vector3d p = Vector3d(sin_alp1 * (this->OR.GV[0].r + dx_out * 0.5)+ this->OR.GV[0].Rx, 0, cos_alp1 * (this->OR.GV[0].r + dx_out * 0.5)) + er * this->OR.GV[0].Rr;
	this->BL.set_y(bl_x_out, bl_v, q, w);

	this->BOP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);

	// θvZΜΚΖδr
	double _fs = 100 * (0.4 * 0.2 + 0.6 * 0.3);
	double bp = this->BL.r - dx_out * 0.5;
	Vector3d _Fbs = Vector3d(0, -_fs, 0);
	Vector3d _Tbs = Vector3d(_fs *bp * cos_alp1, 0, -_fs * bp * sin_alp1);
	double err = 0.05;
	EXPECT_NEAR((Fbs - _Fbs).norm(), 0, _Fbs.norm()*err);
	EXPECT_NEAR((Tbs - _Tbs).norm(), 0, _Tbs.norm()*err);


	return;
}

// (2) ΚͺOΦΙΞ΅ΔXsϋόΙρ]΅½ΜθCπvZ
// ½Ύ΅CN[CWΜέΖ΅½
TEST_F(B4P_BallRingPairTest, calc_Sliding_2) {
	this->init_25BSWZ();

	int i = 0;	// aΤ
	double a = 0.668610936573505e-3; // ΪGΘ~·a[m]
	double Pmax = 1;		// ΕεΚ³igνΘ’j
	double F_norm = 100;	// ]?ΜΧd[N]
	double fratio = 0.0;	// ϋΪG
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double sin_alp1 = 0.540517404112816;		// OΦΪGpΜ³·(balacvZΚ)
	double cos_alp1 = 0.841332832980588;		// OΦΪGpΜ]·(balacvZΚ)
	double dx_out = 4.389357066845456E-006;		// OΦ€ΪίΚ(balacvZΚ)


	this->BOP.TR = new Tribology::Stab_Traction();	// gNVWπ0.2ΙΕθ
	this->BOP.CL = new Tribology::Stab_Coulomb();	// N[CWπ0.3ΙΕθ

	// Κͺ¬x0ΕCΪGΚΌϋόπ²Ιρ]΅Δ’ικπΌθ
	Vector3d x = Vector3d(0, 0, 0);
	Vector3d v = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d w = Vector3d(0, 0, 0);
	this->OR.set_y(x, v, q, w);
	Vector3d er = Vector3d(0, 0, 1);
	double R_out = this->OR.GV[0].r + dx_out - this->BL.r;
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + this->OR.GV[0].Rx, 0, cos_alp1 * R_out) + er * this->OR.GV[0].Rr;
	Vector3d bl_w = Vector3d(sin_alp1, 0, cos_alp1);
	Vector3d p = Vector3d(sin_alp1 * (this->OR.GV[0].r + dx_out * 0.5) + this->OR.GV[0].Rx, 0, cos_alp1 * (this->OR.GV[0].r + dx_out * 0.5)) + er * this->OR.GV[0].Rr;
	this->BL.set_y(bl_x_out, v, q, bl_w);
	this->BOP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);



	// θvZΜΚΖδr
	double _fs = 100 * (0.4 * 0.2 + 0.6 * 0.3);
	double bp = this->BL.r - dx_out * 0.5;
	Vector3d _Fbs = Vector3d(0, 0, 0);
	EXPECT_NEAR((Fbs - _Fbs).norm(), 0, 0.1);			// θCΜg[^ͺ0
	Vector3d eg = Vector3d(sin_alp1, 0, cos_alp1);
	double dir = Tbs.dot(eg);
	Vector3d Dir = Tbs.cross(eg);
	EXPECT_LT(dir, 0);								// gNΜό«ͺΪGΚΙΞ΅Δ{[€ΙΘΑΔ’ι±ΖπmF
	EXPECT_NEAR(Dir.norm(), 0, 0.1);				// gNͺΪGΚΙΞ΅ΔΌΕ ι©mF
	return;
}

// (3) yOΦ+€zΚΪGpͺ30ΕΩΪ]ͺθ΅Δ’ιΖ«ΜCΝΜό«π]Ώ
TEST_F(B4P_BallRingPairTest, calc_Sliding_3) {
	this->init_50BSWZ();
	int i = 1;							// aΤ
	double a = 0.3160272198e-3;			// ΪGΘ~·a[m]
	double Pmax = 2e6;					// ΕεΚ³igνΘ’j
	double F_norm = 20;		// ]?ΜΧd[N]
	double fratio = 0.0;				// ϋΪG
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double sin_alp1 = 0.5;				// 30Μsin
	double cos_alp1 = sqrt(3) * 0.5;	// 30Μcos
	double dx_out = 0;		// OΦ€ΪίΚ(balacvZΚ)

	// ¬ΜFΘ΅CN[CFtanJ[u
	this->BOP.TR = new Tribology::Stab_Traction();
	this->BOP.CL = new Tribology::Tangent();
	this->BOP.FT = new Tribology::FilmThicknessNothing();
	// Κͺ]ͺθ(vρθ)ΕCΪGΚΌϋόπ²Ιρ]΅Δ’ικπΌθ
	double wx = 20;
	Vector3d bl_v = Vector3d(0, wx * (-dx_out + this->BL.r), 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d bl_w = Vector3d(wx * cos_alp1, 0, -wx * sin_alp1);

	// ΚΚuΝaRSΜΚuπξΙθ
	double R_out = this->OR.GV[i].r + dx_out - this->BL.r;
	Vector3d er = Vector3d(0, 0, 1);
	Vector3d Ro = Vector3d(this->OR.GV[i].Rx, 0, this->OR.GV[i].Rr);
	Vector3d bl_x = Ro + Vector3d(sin_alp1 * R_out, 0, cos_alp1 * R_out);
	this->BL.set_y(bl_x, bl_v, q, bl_w);
	Vector3d zero = Vector3d::Zero();
	this->OR.set_y(zero, zero, q, zero);
	double rr = this->OR.GV[i].r + dx_out * 0.5;
	Vector3d p = Ro + Vector3d(sin_alp1 * rr, 0, cos_alp1 * rr);

	this->BOP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);



	// CΝE[gΜ]Ώ
	Vector3d e_fs = Fbs.normalized();
	Vector3d e_ts = Tbs.normalized();
	Vector3d e_fsex(0, -1, 0);
	Vector3d e_tsex(cos_alp1, 0, -sin_alp1);
	EXPECT_NEAR((e_fs - e_fsex).norm(), 0, 1e-6);			// θCΜό«
	EXPECT_NEAR((e_ts - e_tsex).norm(), 0, 1e-6);			// [gΜό«

	return;
}


// (4) yOΦ+€zΚΪGpͺ0ΕΚ©]²ͺx²ξΕ-45πό’Δ’ικΜΚ[gΜό«ͺ³΅­vZΕ«Δ’ι©Ψ
TEST_F(B4P_BallRingPairTest, calc_Sliding_4) {
	this->init_50BSWZ();
	int i = 1;							// aΤ
	double a = 0.3160272198e-3;			// ΪGΘ~·a[m]
	double Pmax = 2e6;					// ΕεΚ³igνΘ’j
	double F_norm = 20;					// ]?ΜΧd[N]
	double fratio = 0.0;				// ϋΪG
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double dx_out = 0;					// OΦ€ΪίΚ
	double cos_45 = 1.0 / sqrt(2.0); 
	double sin_45 = 1.0 / sqrt(2.0);

	this->BOP.TR = new Tribology::Stab_Traction();
	this->BOP.CL = new Tribology::Tangent();
	this->BOP.FT = new Tribology::FilmThicknessNothing();
	// Κͺ¬x0ΕCΪGΚΌϋόπ²Ιρ]΅Δ’ικπΌθ
	double wx = 100;
	Vector3d bl_v = Vector3d(0, wx * (this->BL.r + dx_out * 0.5) * cos_45, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d bl_w = Vector3d(wx * sin_45, 0,  -wx * cos_45);

	// ΚΚuΝaRSΜΚuπξΙθ
	double R_out = this->OR.GV[i].r + dx_out - this->BL.r;
	Vector3d bl_x = Vector3d(this->OR.GV[i].Rx, 0, R_out + this->OR.GV[i].Rr);
	this->BL.set_y(bl_x, bl_v, q, bl_w);
	Vector3d zero = Vector3d::Zero();
	this->OR.set_y(zero, zero, q, zero);
	double rr = this->OR.GV[i].r + dx_out * 0.5;
	Vector3d p = Vector3d(this->OR.GV[i].Rx, 0, rr + this->OR.GV[i].Rr);


	this->BOP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);

	// θvZΜΚΖδr
	EXPECT_GT(Tbs.z(), 0);								// z+ϋόΙ[gͺ­­
	EXPECT_LT(Tis.z(), 0);								// z-ϋόΙ[gͺ­­

	return;
}

// (5) yOΦ+€zΚΪGpͺ0ΕΚ©]²ͺx²ξΕ+45πό’Δ’ικΜΚ[gΜό«ͺ³΅­vZΕ«Δ’ι©Ψ
TEST_F(B4P_BallRingPairTest, calc_Sliding_5) {
	this->init_50BSWZ();
	int i = 1;							// aΤ
	double a = 0.3160272198e-3;			// ΪGΘ~·a[m]
	double Pmax = 2e6;					// ΕεΚ³igνΘ’j
	double F_norm = 20;		// ]?ΜΧd[N]
	double fratio = 0.0;				// ϋΪG
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double dx_out = 0;		// OΦ€ΪίΚ
	double cos_45 = 1.0 / sqrt(2.0);
	double sin_45 = 1.0 / sqrt(2.0);

	this->BOP.TR = new Tribology::Stab_Traction();
	this->BOP.CL = new Tribology::Tangent();
	this->BOP.FT = new Tribology::FilmThicknessNothing();
	// Κͺ¬x0ΕCΪGΚΌϋόπ²Ιρ]΅Δ’ικπΌθ
	double wx = 100;
	Vector3d bl_v = Vector3d(0, wx * (this->BL.r + dx_out * 0.5) * cos_45, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d bl_w = Vector3d(wx * sin_45, 0, wx * cos_45);

	// ΚΚuΝaRSΜΚuπξΙθ
	double R_out = this->OR.GV[i].r + dx_out - this->BL.r;
	Vector3d bl_x = Vector3d(this->OR.GV[i].Rx, 0, R_out + this->OR.GV[i].Rr);
	this->BL.set_y(bl_x, bl_v, q, bl_w);
	Vector3d zero = Vector3d::Zero();
	this->OR.set_y(zero, zero, q, zero);
	double rr = this->OR.GV[i].r + dx_out * 0.5;
	Vector3d p = Vector3d(this->OR.GV[i].Rx, 0, rr + this->OR.GV[i].Rr);


	this->BOP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);
	// θvZΜΚΖδr
	EXPECT_LT(Tbs.z(), 0);								// z-ϋόΙ[gͺ­­

	return;
}

// (6) yΰΦ-€zΚΪGpͺ0ΕΚ©]²ͺx²ξΕ+45πό’Δ’ικΜΚ[gΜό«ͺ³΅­vZΕ«Δ’ι©Ψ
TEST_F(B4P_BallRingPairTest, calc_Sliding_6) {
	this->init_50BSWZ();
	int i = 0;							// aΤ
	double a = 0.3160272198e-3;			// ΪGΘ~·a[m]
	double Pmax = 2e6;					// ΕεΚ³igνΘ’j
	double F_norm = 20;		// ]?ΜΧd[N]
	double fratio = 0.0;				// ϋΪG
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double dxi = 0;		// ΰΦ€ΪίΚ
	double cos_45 = 1.0 / sqrt(2.0);
	double sin_45 = 1.0 / sqrt(2.0);
	double cos_m135 = -1.0 / sqrt(2.0);
	double sin_m135 = -1.0 / sqrt(2.0);
	this->BIP.TR = new Tribology::Stab_Traction();
	this->BIP.CL = new Tribology::Tangent();
	this->BIP.FT = new Tribology::FilmThicknessNothing();
	// ΪG_Ι¨―ιθ¬xͺ0ΜπΌθ
	double wx = 100;
	Vector3d bl_v = Vector3d(0, wx * (this->BL.r + dxi * 0.5) * cos_45, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d bl_w = Vector3d(wx * sin_45, 0, wx * cos_45);

	// ΚΚuΝaRSΜΚuπξΙθ
	double R_in = this->IR.GV[i].r + dxi - this->BL.r;
	Vector3d bl_x = Vector3d(this->IR.GV[i].Rx, 0, this->IR.GV[i].Rr - R_in);
	this->BL.set_y(bl_x, bl_v, q, bl_w);
	Vector3d zero = Vector3d::Zero();
	this->IR.set_y(zero, zero, q, zero);
	double rr = this->IR.GV[i].r + dxi * 0.5;
	Vector3d p = Vector3d(this->IR.GV[i].Rx, 0, rr + this->IR.GV[i].Rr);
	this->BIP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);
	EXPECT_LT(Tbs.z(), 0);								// z-ϋόΙ[gͺ­­

	return;
}


// (1) ³ΜΪGpΜΜOΦΪGpπvZ΅DθvZΜlΖδr
TEST_F(B4P_BallRingPairTest, ContactAngle_1) {
	this->init_25BSWZ();

	double err = 0.01;
	double sin_alp1 = 0.5;				// OΦΪGpΜ³·
	double cos_alp1 = 0.8660254;		// OΦΪGpΜ]·
	double alp = this->BOP.ContactAngle(cos_alp1, sin_alp1);
	double _alp = Unit::deg2rad(30);
	EXPECT_NEAR(_alp, alp, alp * err);
	return;
}

// (2) ΜΪGpΜΜOΦΪGpπvZ΅DθvZΜlΖδr
TEST_F(B4P_BallRingPairTest, ContactAngle_2) {
	this->init_25BSWZ();

	double err = 0.01;
	double sin_alp1 = -0.5;				// OΦΪGpΜ³·
	double cos_alp1 = 0.8660254;		// OΦΪGpΜ]·
	double alp = this->BOP.ContactAngle(cos_alp1, sin_alp1);
	double _alp = Unit::deg2rad(-30);
	EXPECT_NEAR(_alp, alp, abs(alp) * err);
	return;
}

// (3) ³ΜΪGpΜΜΰΦΪGpπvZ΅DθvZΜlΖδr
TEST_F(B4P_BallRingPairTest, ContactAngle_3) {
	this->init_25BSWZ();

	double err = 0.01;
	double sin_alp1 = 0.5;				// OΦΪGpΜ³·
	double cos_alp1 = -0.8660254;		// OΦΪGpΜ]·
	double alp = this->BIP.ContactAngle(cos_alp1, sin_alp1);
	double _alp = Unit::deg2rad(30);
	EXPECT_NEAR(_alp, alp, alp * err);
	return;
}

// (4) ΜΪGpΜΜΰΦΪGpπvZ΅DθvZΜlΖδr
TEST_F(B4P_BallRingPairTest, ContactAngle_4) {
	this->init_25BSWZ();

	double err = 0.01;
	double sin_alp1 = -0.5;				// OΦΪGpΜ³·
	double cos_alp1 = -0.8660254;		// OΦΪGpΜ]·
	double alp = this->BIP.ContactAngle(cos_alp1, sin_alp1);
	double _alp = Unit::deg2rad(-30);
	EXPECT_NEAR(_alp, alp, abs(alp) * err);
	double alp__ = Unit::rad2deg(atan2(sin_alp1, cos_alp1));
	return;
}

// (1) πͺςΜmF
TEST_F(B4P_BallRingPairTest, init_Tribology_1) {
	this->init_25BSWZ();

	FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::Aihara;
	this->BOP.init_Tribology(FI.TB);
	Tribology::AiharaR* aih = dynamic_cast<Tribology::AiharaR*>(this->BOP.RR);
	EXPECT_TRUE(aih != nullptr);

	FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::Fujiwara;
	this->BOP.init_Tribology(FI.TB);
	Tribology::Fujiwara* fj = dynamic_cast<Tribology::Fujiwara*>(this->BOP.RR);
	EXPECT_TRUE(fj != nullptr);

	FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::Houpert;
	this->BOP.init_Tribology(FI.TB);
	Tribology::Houpert* hou = dynamic_cast<Tribology::Houpert*>(this->BOP.RR);
	EXPECT_TRUE(hou != nullptr);

	FI.TB.coulomb = B4P_In::Tribology::Coulomb::CoulombNothing;
	this->BOP.init_Tribology(FI.TB);
	Tribology::CoulombNothing* cln = dynamic_cast<Tribology::CoulombNothing*>(this->BOP.CL);
	EXPECT_TRUE(cln != nullptr);

	FI.TB.coulomb = B4P_In::Tribology::Coulomb::Tangent;
	this->BOP.init_Tribology(FI.TB);
	Tribology::Tangent* tan = dynamic_cast<Tribology::Tangent*>(this->BOP.CL);
	EXPECT_TRUE(tan != nullptr);

	FI.TB.filmThickness = B4P_In::Tribology::FilmThickness::FilmThicknessNothing;
	this->BOP.init_Tribology(FI.TB);
	Tribology::FilmThicknessNothing* flt = dynamic_cast<Tribology::FilmThicknessNothing*>(this->BOP.FT);
	EXPECT_TRUE(flt != nullptr);

	FI.TB.filmThickness = B4P_In::Tribology::FilmThickness::HamrockDowsonHc;
	this->BOP.init_Tribology(FI.TB);
	Tribology::HamrockDowsonHc* hmd = dynamic_cast<Tribology::HamrockDowsonHc*>(this->BOP.FT);
	EXPECT_TRUE(hmd != nullptr);
	return;
}

