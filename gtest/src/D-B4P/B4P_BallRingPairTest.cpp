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

	// 軸受緒元を50BSWZ02のものに設定
	void init_50BSWZ() {
		// 潤滑
		FI.LB.eta0 = 0.5;	// 粘度（Pa*s）(適当)
		FI.LB.beta0 = 0.5;	// 温度粘度係数（適当）
		FI.LB.k0 = 0.145;	// 油熱伝導率[W/(m・K)]
		FI.LB.alpha0 = 0.2;	// 圧力粘度係数[mm2/kgf]
		FI.LB.lm0 = -1;	// メニスカス長さ


		// 玉物性値
		double D = 0.0111;				// 玉径[m]
		double E = 208000000000;		// ヤング率[Pa]
		double por = 0.29;				// ポアソン比[-]
		double den = 7830;				// 密度[kg/m^3]
		double rms = 0.00000002;		// 粗さrms[m]
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}

		this->BL.init(D, E, por, den, rms, x_const, Rx_const);
		int msmax = 21;
		FO.allocate(1, _MAX_CONTACT_, msmax);

		// 内輪物性値
		FI.ballpcd = 0.0675;
		double clr = 0.000008;
		FI.IR.R[0] = 0.00589;
		FI.IR.R[1] = 0.00589;
		FI.IR.Rox[0] = 0.000167;
		FI.IR.Rox[1] = -0.000167;
		// 溝R中心PCD[m]
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
		// 潤滑理論
		FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::RollingResistanceNothing;
		FI.TB.coulomb = B4P_In::Tribology::Tangent;
		FI.TB.coulomb_slope = 1e100;
		FI.TB.hysteresis = B4P_In::Tribology::HysteresisNothing;
		this->IR.init(FI.IR, FI.ballpcd);
		this->BIP.link(&this->BL, &this->IR);
		this->BIP.init(FI);

		// 外輪物性値
		FI.OR.Rox[0] = 0.0001665;
		FI.OR.Rox[1] = -0.0001665;
		FI.OR.R[0] = 0.00589;
		FI.OR.R[1] = 0.00589;
		// 溝R中心PCD[m]
		FI.OR.Rod[0] = FI.ballpcd + clr * 0.5 - 2 * sqrt((FI.OR.R[0] - D * 0.5)*(FI.OR.R[0] - D * 0.5) - FI.OR.Rox[0] * FI.OR.Rox[0]);
		FI.OR.Rod[1] = FI.ballpcd + clr * 0.5 - 2 * sqrt((FI.OR.R[1] - D * 0.5)*(FI.OR.R[1] - D * 0.5) - FI.OR.Rox[1] * FI.OR.Rox[1]);
		FI.OR.hedge[0];
		FI.OR.hedge[1];
		FI.OR.rms = 0.00000006;
		FI.OR.m = 1;
		FI.OR.E = 207760000000;
		FI.OR.por = 0.29;
		FI.OR.Ix, FI.OR.Iyz, FI.OR.Iyz;
		FI.LB.eta0 = 0.5;	// 粘度（Pa*s）(適当)
		FI.LB.beta0 = 0.5;	// 温度粘度係数（適当）
		FI.LB.k0 = 0.145;	// 油熱伝導率[W/(m・K)]
		FI.LB.alpha0 = 0.2;	// 圧力粘度係数[mm2/kgf]
		FI.LB.lm0 = -1;	// メニスカス長さ
		FI.BOP.mu = 0.1;	// 摩擦係数
		FI.BOP.dzeta = 0.2;	// 減衰係数
		this->OR.init(FI.OR, FI.ballpcd);
		this->BOP.link(&this->BL, &this->OR);
		this->BOP.init(FI);


		return;
	};

	// 軸受緒元を25BSWZ01のものに設定
	void init_25BSWZ() {
		// 入力データは下記のものを用いた
		// \\ans00978\kiken\BRAIN\ソースコード他\DBRAINpp\docs\06_D4Bv100\05_単体テスト\B4P_BallRingPair
		// 軸受型番　25BSWZ01
		// Pairクラス物性値
		this->FI.BOP.dzeta = 0.2;	// 減衰項
		this->FI.BOP.mu = 0.1;		// 摩擦係数
		this->FI.BIP.dzeta = 0.2;	// 減衰項
		this->FI.BIP.mu = 0.1;		// 摩擦係数

		// 玉物性値
		double D = 0.00635;				// 玉径[m]
		double E = 207760000000;		// ヤング率[Pa]
		double por = 0.29;				// ポアソン比[-]
		double den = 7830;				// 密度[kg/m^3]
		double rms = 0.00000002;		// 粗さrms[m]
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}
		this->BL.init(D, E, por, den, rms, x_const, Rx_const);
		int msmax = 21;
		FO.allocate(1, _MAX_CONTACT_, msmax);

		// 内輪物性値
		FI.ballpcd = 0.0355;
		double clr = 0;
		FI.IR.Rox[0] = 0.0000725;
		FI.IR.Rox[1] = -0.0000725;
		FI.IR.R[0] = 0.00332105;
		FI.IR.R[1] = 0.00332105;
		// 溝R中心PCD[m]（溝pcd半径，0.0357535 m くらい）
		FI.IR.Rod[0] = FI.ballpcd - clr * 0.5 + 2 * sqrt((FI.IR.R[0] - D * 0.5)*(FI.IR.R[0] - D * 0.5) - FI.IR.Rox[0] * FI.IR.Rox[0]);
		FI.IR.Rod[1] = FI.ballpcd - clr * 0.5 + 2 * sqrt((FI.IR.R[1] - D * 0.5)*(FI.IR.R[1] - D * 0.5) - FI.IR.Rox[1] * FI.IR.Rox[1]);
		FI.IR.hedge[0];
		FI.IR.hedge[1];
		FI.IR.rms = 0.00000006;
		FI.IR.m = 1;
		FI.IR.E = 207760000000;
		FI.IR.por = 0.29;
		FI.IR.Ix, FI.IR.Iyz, FI.IR.Iyz;

		// 潤滑理論
		FI.TB.rollingresistance = B4P_In::Tribology::RollingResistance::RollingResistanceNothing;
		FI.TB.coulomb = B4P_In::Tribology::Tangent;
		FI.TB.coulomb_slope = 1000;
		FI.TB.hysteresis = B4P_In::Tribology::HysteresisNothing;
		this->IR.init(FI.IR, FI.ballpcd);
		this->BIP.link(&this->BL, &this->IR);
		this->BIP.init(FI);

		// 外輪物性値
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
		FI.LB.eta0 = 0.5;	// 粘度（Pa*s）(適当)
		FI.LB.beta0 = 0.5;	// 温度粘度係数（適当）
		FI.LB.k0 = 0.145;	// 油熱伝導率[W/(m・K)]
		FI.LB.alpha0 = 0.2;	// 圧力粘度係数[mm2/kgf]
		FI.LB.lm0 = -1;	// メニスカス長さ
		this->OR.init(FI.OR, FI.ballpcd);
		this->BOP.link(&this->BL, &this->OR);
		this->BOP.init(FI);

	};

	virtual void SetUp() {
		Vector3d x = Vector3d(0, 0, 0);
		Vector3d v = Vector3d(0, 0, 0);
		Quaterniond q = Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）
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



// (1) 保持器・玉・内外輪自転速度をパラメータとした式を作成し，比較．
TEST_F(B4P_BallRingPairTest, get_us_ur_1) {
	this->init_25BSWZ();
	double wc = 44.4604328912874;		// 保持器公転速度[rad/s]
	double wrn_0 = 0.0;					// 外輪回転速度[rad/s]
	double wrn_1 = 104.719755119660;	// 内輪回転速度[rad/s]
	double omega_x = -223.906546549487;	// 玉自転速度 ωx[rad/s]
	double omega_z = 287.884790174018;  // 玉自転速度 ωz[rad/s]
	double dx_out = 4.389357066845456E-006;		// 外輪側接近量(balac計算結果)
	double dx_in = 4.518845377961239E-006;		// 内輪側接近量(balac計算結果)
	double sin_alp1 = 0.537033675917149;		// 外輪接触角の正弦(balac計算結果)
	double sin_alp2 = 0.537125019737040;		// 内輪接触角の正弦(balac計算結果)
	double cos_alp1 = 0.843560804525029;		// 外輪接触角の余弦(balac計算結果)
	double cos_alp2 = -0.843502645622694;		// 内輪接触角の余弦(balac計算結果)

	// 玉相対位置は接近量から逆算
	double R_out = this->OR.GV[0].r - this->BL.r;
	double R_in = this->IR.GV[0].r - this->BL.r;
	Vector3d er = Vector3d(0, 0, 1);			// 径方向ベクトル
	Vector3d bl_x_in = Vector3d(0, 0, 0);
	Vector3d bl_x_out = Vector3d(0, 0, 0);

	Vector3d x = Vector3d(0, 0, 0);
	Vector3d v = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）
	Vector3d wo = Vector3d(wrn_0, 0, 0);
	this->OR.set_y(x, v, q, wo);
	Vector3d wi = Vector3d(wrn_1, 0, 0);
	this->IR.set_y(x, v, q, wi);
	Vector3d bl_x = Vector3d(0, 0, 0.5 * FI.ballpcd);
	Vector3d bl_v = Vector3d(0, -0.5 * FI.ballpcd * wc, 0);
	Vector3d bl_w = Vector3d(omega_x, 0, omega_z);
	this->BL.set_y(bl_x, bl_v, q, bl_w);

	// 玉中心に対する接触点の相対距離
	Vector3d _p_out = Vector3d(this->BL.r * sin_alp1, 0, this->BL.r * cos_alp1 + 0.5 * FI.ballpcd);
	Vector3d _p_in = Vector3d(this->BL.r * sin_alp2, 0, this->BL.r * cos_alp2 + 0.5 * FI.ballpcd);
	Vector3d uso_, uro_, usi_, uri_;
	this->BOP.get_us_ur(_p_out, uso_, uro_);

	double uso_norm = uso_.norm(), uro_norm = uro_.norm();
	this->BIP.get_us_ur(_p_in, usi_, uri_);
	double usi_norm = usi_.norm(), uri_norm = uri_.norm();

	// baltacに実装されていた式を参考に自作した式(一部変数の定義方法が異なるため改変)
	// （接触点の速度の計算方法も本プログラムの計算方法に合わせこんでいる）
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
	// 許容誤差
	EXPECT_NEAR(uro_norm, abs(baltac_Uo), abs(baltac_Uo)*err);
	EXPECT_NEAR(uso_norm, abs(baltac_uso), abs(baltac_uso)*err);
	EXPECT_NEAR(uri_norm, abs(baltac_Ui), abs(baltac_Ui)*err);
	EXPECT_NEAR(usi_norm, abs(baltac_usi), abs(baltac_usi)*err);

	return;
}

// (1) 保持器・玉・内外輪自転速度をパラメータとした式を作成し，比較．
TEST_F(B4P_BallRingPairTest, get_us_ur2_1) {

	this->init_25BSWZ();
	double wc = 44.4604328912874;		// 保持器公転速度[rad/s]
	double wrn_0 = 0.0;					// 外輪回転速度[rad/s]
	double wrn_1 = 104.719755119660;	// 内輪回転速度[rad/s]

	double omega_x = -223.906546549487;	// 玉自転速度 ωx[rad/s]
	double omega_z = 287.884790174018;  // 玉自転速度 ωz[rad/s]
	double dx_out = 4.389357066845456E-006;		// 外輪側接近量(balac計算結果)
	double dx_in = 4.518845377961239E-006;		// 内輪側接近量(balac計算結果)
	double sin_alp1 = 0.537033675917149;		// 外輪接触角の正弦(balac計算結果)
	double sin_alp2 = 0.537125019737040;		// 内輪接触角の正弦(balac計算結果)
	double cos_alp1 = 0.843560804525029;		// 外輪接触角の余弦(balac計算結果)
	double cos_alp2 = -0.843502645622694;		// 内輪接触角の余弦(balac計算結果)

	// 玉相対位置は接近量から逆算
	double R_out = this->OR.GV[0].r - this->BL.r;
	double R_in = this->IR.GV[0].r - this->BL.r;
	Vector3d er = Vector3d(0, 0, 1);			// 径方向ベクトル
	Vector3d bl_x_in = Vector3d(0, 0, 0.5 * FI.ballpcd);
	Vector3d bl_x_out = Vector3d(0, 0, 0.5 * FI.ballpcd);

	Vector3d x = Vector3d(0, 0, 0);
	Vector3d v = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0); // 初期クォータニオンは全て [0,0,0,1] とする．（ただしコンストラクタの仕様上順番が異なる）
	Vector3d wo = Vector3d(wrn_0, 0, 0);
	this->OR.set_y(x, v, q, wo);
	Vector3d wi = Vector3d(wrn_1, 0, 0);
	this->IR.set_y(x, v, q, wi);
	Vector3d bl_v = Vector3d(0, -0.5 * FI.ballpcd * wc, 0);
	Vector3d bl_w = Vector3d(omega_x, 0, omega_z);
	this->BL.set_y(bl_x_in, bl_v, q, bl_w);

	// 玉中心に対する接触点の相対距離
	Vector3d _p_out = Vector3d(this->BL.r * sin_alp1, 0, this->BL.r * cos_alp1 + 0.5 * FI.ballpcd);
	Vector3d _p_in = Vector3d(this->BL.r * sin_alp2, 0, this->BL.r * cos_alp2 + 0.5 * FI.ballpcd);
	Vector3d uso_, uro_, usi_, uri_;
	this->BOP.get_us_ur2(_p_out, er, uso_, uro_);

	double uso_norm = uso_.norm(), uro_norm = uro_.norm();
	this->BIP.get_us_ur2(_p_in, er, usi_, uri_);
	double usi_norm = usi_.norm(), uri_norm = uri_.norm();

	// baltacに実装されていた式を参考に自作した式(一部変数の定義方法が異なるため改変)
	// （接触点の速度の計算方法も本プログラムの計算方法に合わせこんでいる）
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
	// 許容誤差
	EXPECT_NEAR(uro_norm, abs(baltac_Uo), abs(baltac_Uo)*err);
	EXPECT_NEAR(uso_norm, abs(baltac_uso), abs(baltac_uso)*err);
	EXPECT_NEAR(uri_norm, abs(baltac_Ui), abs(baltac_Ui)*err);
	EXPECT_NEAR(usi_norm, abs(baltac_usi), abs(baltac_usi)*err);

	return;
}


// (1) 摩擦係数を固定値にしたときに玉ー内外輪間の荷重が手計算と同じになっているか検証
TEST_F(B4P_BallRingPairTest, calc_force_1) {
	this->init_25BSWZ();

	// 接触状態
	double sin_alp1 = 0.540517404112816;		// 外輪接触角の正弦(balac計算結果)
	double cos_alp1 = 0.841332832980588;		// 外輪接触角の余弦(balac計算結果)
	double dx_out = 4.389357066845456E-006;		// 外輪側接近量(balac計算結果)
	double _Qo = 21.3811947318042;				// 転動体荷重[kgf](balac計算結果)

	this->BOP.TR = new Tribology::Stab_Traction();	// トラクション係数を0.2に固定
	this->BOP.CL = new Tribology::Stab_Coulomb();	// クーロン摩擦係数を0.3に固定
	this->BOP.RR = new Tribology::Stab_RollingResist();	// 転がり粘性抵抗を定数に固定
	this->BOP.FT = new Tribology::FilmThicknessNothing();	// 油膜厚さを0に固定

	// 玉が無回転で並進運動している場合を仮定
	int i = 1;			// 接触する溝番号
	Vector3d x = Vector3d(0, 0, 0), v = Vector3d(0, 0, 0), w = Vector3d(0, 0, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	this->OR.set_y(x, v, q, w);
	Vector3d er = Vector3d(0, 0, 1);
	double dx = this->OR.GV[i].Rx;
	double R_out = this->OR.GV[i].r + dx_out - this->BL.r;			// 溝R中心から玉中心の距離
	double rp = this->OR.GV[i].r + dx_out * 0.5;					// 溝R中心から接触点の距離
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + dx, 0, cos_alp1 * R_out) + er * this->OR.GV[i].Rr;
	Vector3d bl_v = Vector3d(0, -1, 0);
	this->BL.set_y(bl_x_out, bl_v, q, w);


	Vector3d Fbi, Tbi, Fib, Tib;

	this->BOP.calc_force(Fbi, Tbi, Fib, Tib);

	// 手計算の結果と比較
	double _fn = _Qo * 9.8;
	Vector3d _Fn = Vector3d(-sin_alp1, 0, -cos_alp1) * _fn;	// 垂直荷重[N]
	Vector3d _Fs = Vector3d(0, _fn * 0.3, 0);				// 滑り摩擦[N]
	Vector3d _Fbi = _Fn + _Fs;
	double err = 0.05;
	EXPECT_NEAR((Fbi - _Fbi).norm(), 0, _Fbi.norm() * err);
	double _rb = this->BL.r - dx_out * 0.5;
	double _Ts_norm = _fn * 0.3 * _rb;
	Vector3d _Ts = Vector3d(-cos_alp1, 0, sin_alp1) * _Ts_norm;			// 滑り摩擦[N]
	Vector3d _Tr = Vector3d(cos_alp1, 0, -sin_alp1) * 0.1;
	Vector3d _Tbi = _Ts + _Tr;
	EXPECT_NEAR((Tbi - _Tbi).norm(), 0, _Tbi.norm() * err);

	// write()のテストを兼ねて接触楕円や接触角を検証

	this->BOP.write(FO.BOP[0]);
	Vector3d Fn(FO.BOP[0].GV[i].Fn);		// 垂直荷重[N]
	EXPECT_NEAR((Fn - _Fn).norm(), 0, _Fn.norm() * err);
	Vector3d Fs(FO.BOP[0].GV[i].Fs);		// 滑り摩擦[N]
	EXPECT_NEAR((Fs - _Fs).norm(), 0, _Fs.norm() * err);
	double ea = FO.BOP[0].GV[i].a;											// 楕円長半径[mm](実際)
	double _ea = 0.668610936573505e-3;										// 楕円長半径[mm](期待値)
	EXPECT_NEAR(_ea, ea, _ea * err);
	double eb = FO.BOP[0].GV[i].b;											// 楕円短半径[mm](実際)
	double _eb = 9.720461846754701e-5;								// 楕円短半径[mm](期待値)
	EXPECT_NEAR(_eb, eb, _eb * err);
	double alp = FO.BOP[0].GV[i].phi;											// 接触角[°](実際)
	double _alp = Unit::deg2rad(32.7188677338554);										// 接触角[°](期待値)
	EXPECT_NEAR(_alp, alp, _alp * err);

	return;
}

// (1) 接触角0°における玉接近量と各種方向ベクトルを確認
TEST_F(B4P_BallRingPairTest, how_Contact_1) {
	this->init_25BSWZ();

	// 接触角0°・接近量1mmのときの玉中心位置を逆算
	// 玉中心y座標 = 接近量 + 溝R半径 - 玉半径 + 溝R中心y座標
	// ただし，溝中心位置は(0.0000725, 0.035246430660370774* 0.5, 0)
	// 接触する溝は +X 側とする

	Vector3d x(-0.0000725, (0.00332105 + 0.035246430660370774 * 0.5 - 0.00635 * 0.5 + 0.001), 0);

	Vector3d er, eg;
	double dx;

	int i = 1; 
	bool c = this->BOP.how_Contact(i, x, er, eg, dx);

	double err = 0.01;					// 許容誤差
	Vector3d expcted_er(0, 1, 0);		// 玉位相ベクトル
	Vector3d expcted_eg(0, 1, 0);		// 玉接触方向ベクトル
	double expected_dx = 0.001;			// 玉接近量

	EXPECT_NEAR((er - expcted_er).norm(), 0, er.norm()*err);
	EXPECT_NEAR((eg - expcted_eg).norm(), 0, eg.norm()*err);
	EXPECT_NEAR(dx, expected_dx, dx*err);

	return;
};

// (2) ボール中心がトーラスから出ていたら接触なしの判定となることを確認
TEST_F(B4P_BallRingPairTest, how_Contact_2) {
	this->init_25BSWZ();

	// ボール中心をトーラスより外側に設定，接触する溝は+X側に設定
	Vector3d x(-0.0000725, 1.0, 0);
	int i = 1;

	Vector3d er, eg;
	double dx;
	bool c = this->BOP.how_Contact(i, x, er, eg, dx);

	EXPECT_EQ(false, c);
	return;
};

// (3) 食い込み量が負値の場合接触なしの判定となることを確認
TEST_F(B4P_BallRingPairTest, how_Contact_3) {
	this->init_25BSWZ();

	// 接近量-10umかつ接触角0°となるように玉位置を指定
	Vector3d x(-0.0000725, (0.00332105 + 0.035246430660370774 * 0.5 - 0.00635 * 0.5 - 0.0001), 0);
	//Vector3d v(0.0, -0.01, 0.0);
	int i = 1;

	Vector3d er, eg;
	double dx;

	bool c = this->BOP.how_Contact(i, x,  er, eg, dx);
	EXPECT_EQ(false, c);
	return;
};


// (4) 玉が溝中心点と重なるとき，接触なしになるか確認
TEST_F(B4P_BallRingPairTest, how_Contact_4) {
	this->init_25BSWZ();

	// 玉位置を溝中心位置に配置
	Vector3d x(-0.0000725, 0.035246430660370774 * 0.5, 0);
	int i = 1;

	Vector3d er, eg;
	double dx;

	// 通常はありえないが，玉径を1.2倍にして玉が溝中心にある時でも外輪と接触するように調整
	double D = 0.00635 * 1.2;		// 玉径[m]
	double E = 207760000000;		// ヤング率[Pa]
	double por = 0.29;				// ポアソン比[-]
	double den = 7830;				// 密度[kg/m^3]
	double rms = 0.00000002;		// 粗さrms[m]
	bool x_const[3], Rx_const[3];
	for (int i = 0; i < 3; i++) {
		x_const[i] = false;
		Rx_const[i] = false;
	}
	this->BL.init(D, E, por, den, rms, x_const, Rx_const);


	bool c = this->BOP.how_Contact(i, x, er, eg, dx);

	double err = 0.05; // 許容誤差
	Vector3d expcted_er(0, 1, 0);
	Vector3d expcted_eg(0, 1, 0);
	double expected_dx = 0.001;

	EXPECT_EQ(false, c);
	return;
};

// (5) 位相角と接触角が非ゼロかつ，接触・非接触の境界付近の挙動を確認
TEST_F(B4P_BallRingPairTest, how_Contact_5) {
	this->init_25BSWZ();

	// 玉を位相角60°，接触角30°，接近量0の位置に配置
	// x座標 = 溝R中心x座標 + (接近量 + 溝R半径 - 玉半径) * sin30°
	// y座標 = ((接近量 + 溝R半径 - 玉半径) * cos30° + 溝PCD半径) * sin60°
	// z座標 = ((接近量 + 溝R半径 - 玉半径) * cos30° + 溝PCD半径) * cos60°

	Vector3d x(5.25E-07, 0.01537169, 0.008874849);
	int i = 1;
	Vector3d er, eg;
	double dx;

	// 接近量0の位置から +x, +y, +z 方向に微小量移動させると接触するか判定
	Vector3d ep_x(1e-6, 0, 0), ep_y(0, 1e-6, 0), ep_z(0, 0, 1e-6);
	bool cx1 = this->BOP.how_Contact(i, x + ep_x, er, eg, dx);
	bool cy1 = this->BOP.how_Contact(i, x + ep_y, er, eg, dx);
	bool cz1 = this->BOP.how_Contact(i, x + ep_z, er, eg, dx);
	EXPECT_EQ(true, cx1);
	EXPECT_EQ(true, cy1);
	EXPECT_EQ(true, cz1);

	// 接近量0の位置から -x, -y, -z 方向に微小量移動させると非接触になるか判定
	bool cx2 = this->BOP.how_Contact(i, x - ep_x, er, eg, dx);
	bool cy2 = this->BOP.how_Contact(i, x - ep_y, er, eg, dx);
	bool cz2 = this->BOP.how_Contact(i, x - ep_z, er, eg, dx);
	EXPECT_EQ(false, cx2);
	EXPECT_EQ(false, cy2);
	EXPECT_EQ(false, cz2);
	return;
};


// (1) 玉速度を変えて計算を行い，減衰を考慮しない場合との比較を行う．
TEST_F(B4P_BallRingPairTest, calc_DynamicHertz_1) {
	this->init_25BSWZ();


	/*****d4b入力条件(軸受緒元はSetUp関数で設定)*****/
	// 玉が軸方向から見て +y 軸方向にある場合を想定．
	// 玉速度が +y 方向（外輪外側に進行）

	int i = 0;									// 溝番号
	double sin_alp1 = 0.540517404112816;		// 外輪接触角の正弦(balac計算結果)
	double cos_alp1 = 0.841332832980588;		// 外輪接触角の余弦(balac計算結果)
	Vector3d er= Vector3d(0, 1, 0);					// 径方向ベクトル
	Vector3d eg_out= Vector3d(sin_alp1, cos_alp1, 0);	// 溝中心→ボール方向ベクトル(baltac接触角から逆算)
	double dx_out= 4.389357066845456E-006;		// 外輪側接近量(balac計算結果)
	double Rx_out, Ry_out, cos_alp, sin_alp, a_out, b_out;
	Vector3d p_out;
	// 玉相対位置は接近量から逆算
	// 玉を位相角 90° (+y方向)，接触角 30°，接近量 4.389 um の位置に配置
	// x座標 = 溝R中心x座標 + (接近量 + 溝R半径 - 玉半径) * sin30°
	// y座標 = ((接近量 + 溝R半径 - 玉半径) * cos30° + 溝PCD半径) * sin90°
	// z座標 = ((接近量 + 溝R半径 - 玉半径) * cos30° + 溝PCD半径) * cos90°
	double R_out = this->OR.GV[0].r + dx_out - this->BL.r;
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + this->OR.GV[0].Rx, cos_alp1 * R_out, 0) + er * this->OR.GV[0].Rr;
	
	// 比較対象となる減衰なしの場合を先に計算．
	double k_out, k_in;
	this->BOP.how_Contact(i, bl_x_out, er, eg_out, dx_out);
	this->BOP.calc_Hertz(i, bl_x_out, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out, k_out);
	double _Qout = k_out * pow(dx_out, 1.5);
	// (1) 玉が接触面に近づく向きに進行する場合，接触力が大きくなるか確認
	Vector3d bl_v1 = Vector3d(1, 0, 0);		// 玉速度（壁面から近づく向き）
	double Qout1 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v1, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_GE(Qout1, _Qout);	

	// (2) 玉が接触面から離れる向きに進行する場合，接触力が小さくなるか確認
	Vector3d bl_v2 = Vector3d(-1, 0, 0);	// 玉速度（壁面から離れる向き）
	double Qout2 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v2, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_LE(Qout2, _Qout);

	// (3) 玉が接触面と並行な向きに進行する場合，接触力が変わらないことを確認(z成分のみ)
	Vector3d bl_v3 = Vector3d(0, 0, -1);	// 玉速度（壁面並行方向）
	double Qout3 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v3, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_NEAR(Qout3, _Qout, 1e-6);

	// (4) 玉が接触面と並行な向きに進行する場合，接触力が変わらないことを確認(xy成分のみ)
	Vector3d bl_v4 = Vector3d(0.841332832980588, -0.540517404112816, 0);	// 玉速度（壁面並行方向）
	double Qout4 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v4, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_NEAR(Qout4, _Qout, 1e-6);

	// (5) 減衰力とHertz接触力の合計値が負の時，0に補正していることを確認
	Vector3d bl_v5 = Vector3d(-1e9, 0, 0);		// 玉速度（壁面から離れる向き）
	double Qout5 = this->BOP.calc_DynamicHertz(i, bl_x_out, bl_v5, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out);
	EXPECT_NEAR(Qout5, 0, 1e-6);	

	return;
};




// (1) 減衰成分がないときのヘルツ接触力を計算し，baltacの計算結果と比較
// baltacと同等の接近量・接触角を入力した条件で，転動体荷重・楕円半径を評価
// （入力条件はアキシアル荷重の場合を用いる．減衰は考慮しない．）
TEST_F(B4P_BallRingPairTest, calc_Hertz_1) {
	this->init_25BSWZ();

	// 軸受緒元はSetUp関数で設定
	// d4b入力条件
	int i = 0;									// 溝番号
	double sin_alp1 = 0.540517404112816;		// 外輪接触角の正弦(balac計算結果)
	double sin_alp2 = 0.540517404112816;		// 内輪接触角の正弦(balac計算結果)
	double cos_alp1 = 0.841332832980588;		// 外輪接触角の余弦(balac計算結果)
	double cos_alp2 = -0.841332832980588;		// 内輪接触角の余弦(balac計算結果)
	Vector3d er = Vector3d(0, 1, 0);			// 径方向ベクトル
	Vector3d eg_out = Vector3d(sin_alp1, cos_alp1, 0); // 溝中心→ボール方向ベクトル(baltac接触角から逆算)
	Vector3d eg_in = Vector3d(sin_alp2, cos_alp2, 0); // 溝中心→ボール方向ベクトル(baltac接触角から逆算)
	double dx_out = 4.389357066845456E-006;		// 外輪側接近量(balac計算結果)
	double dx_in = 4.518845377961239E-006;		// 内輪側接近量(balac計算結果)
	double Rx_out, Ry_out, Rx_in, Ry_in, cos_alp, sin_alp, a_out, b_out, k_out, a_in, b_in, k_in;
	Vector3d p_out, p_in;
	// 玉相対位置は接近量から逆算
	double R_out = this->OR.GV[0].r + dx_out - this->BL.r;
	double R_in = this->IR.GV[0].r + dx_in - this->BL.r;
	Vector3d bl_x_out = Vector3d(sin_alp1 * R_out + this->OR.GV[0].Rx, cos_alp1 * R_out, 0) + er * this->OR.GV[0].Rr;
	Vector3d bl_x_in = Vector3d(sin_alp2 * R_in + this->IR.GV[0].Rx, cos_alp2 * R_in, 0) + er * this->IR.GV[0].Rr;
	// テスト対象関数の呼び出し
	this->BOP.calc_Hertz(i, bl_x_out, er, eg_out, dx_out, Rx_out, Ry_out, p_out, cos_alp, sin_alp, a_out, b_out, k_out);
	double Qout = k_out * pow(dx_out, 1.5);
	this->BIP.calc_Hertz(i, bl_x_in, er, eg_in, dx_in, Rx_in, Ry_in, p_in, cos_alp, sin_alp, a_in, b_in, k_in);
	double Qin = k_in * pow(dx_in, 1.5);

	// 比較対象となるbalac計算結果（デバッグモードにて取得）
	// 1: 外輪，2:内輪
	double ea1 = 0.668610936573505;				// 楕円長半径[mm]
	double ea2 = 0.685074685553727;				// 楕円長半径[mm]
	double eb1 = 9.720461846754701E-002;		// 楕円短半径[mm]
	double eb2 = 8.272063984522840E-002;		// 楕円短半径[mm]
	double _Q1 = 21.3811947318042;				// 転動体荷重[kgf]
	double _Q2 = 21.3811947318041;				// 転動体荷重[kgf]
	double _rho_0 = 2 / 6.35;					// 転動体曲率[1/mm]
	double _rho3_in = -0.301109588834856;		// 溝曲率[1/mm]
	double _rho2_in = 5.579585936574204E-002;	// 軌道の曲率(公転軌道接触半径の逆数）[1/mm]
	double _Rx_in = 1 / (_rho_0 + _rho2_in);
	double _Ry_in = 1 / (_rho_0 + _rho3_in);

	double _rho3_out = -0.301109588834856;
	double _rho2_out = -4.119892685701449E-002;
	double _Rx_out = 1 / (_rho_0 + _rho2_out);
	double _Ry_out = 1 / (_rho_0 + _rho3_out);
	// 接触位置(テストケース)は接近量と接触角から計算
	Vector3d _p_out = Vector3d(sin_alp1 * (this->OR.GV[0].r + dx_out * 0.5) + this->OR.GV[0].Rx, cos_alp1 * (this->OR.GV[0].r + dx_out * 0.5), 0)
		+ er * this->OR.GV[0].Rr;
	Vector3d _p_in = Vector3d(sin_alp2 * (this->IR.GV[0].r + dx_in * 0.5) + this->IR.GV[0].Rx, cos_alp2 * (this->IR.GV[0].r + dx_in * 0.5), 0)
		+ er * this->IR.GV[0].Rr;


	cout << "★外輪側結果" << endl;
	cout << "転動体荷重：" << Qout << "\tbaltac：" << _Q1 * 9.8 << endl;
	cout << "接触楕円長半径：" << a_out << "\tbaltac：" << ea1 / 1000 << endl;
	cout << "接触楕円短半径：" << b_out << "\tbaltac：" << eb1 / 1000 << endl;
	cout << "曲率Rx：" << Rx_out << "\tbaltac：" << _Rx_out / 1000 << endl;
	cout << "曲率Ry：" << Ry_out << "\tbaltac：" << _Ry_out / 1000 << endl;
	cout << "接触点位置誤差：" << (p_out - _p_out).norm() << endl;
	//cout << "p = " << p_out << "\thand：" << _p_out <<endl;
	cout << "★内輪側結果" << endl;
	cout << "転動体荷重：" << Qin << "\tbaltac：" << _Q2 * 9.8 << endl;
	cout << "接触楕円長半径：" << a_in << "\tbaltac：" << ea2 / 1000 << endl;
	cout << "接触楕円短半径：" << b_in << "\tbaltac：" << eb2 / 1000 << endl;
	cout << "曲率Rx：" << Rx_in << "\tbaltac：" << _Rx_in / 1000 << endl;
	cout << "曲率Ry：" << Ry_in << "\tbaltac：" << _Ry_in / 1000 << endl;
	cout << "接触点位置誤差：" << (p_in - _p_in).norm() << endl;
	//cout << "p = " << p_in << "\thand：" << _p_in << endl;
	// baltacではBrew-Hamrockの式ではなく，
	// 無限級数展開による近似で求めているため完全に一致しない．
	// （baltac内部のK(k')とE(k')の計算は運動力学Ⅱp12に掲載）
	// 誤差が2%以内であれば可とした．
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

// (1)玉が内外輪円周方向に移動したときの滑り摩擦を計算し，手計算の数値と比較
// ただし，トラクション係数・クーロン摩擦係数は定数値に固定
TEST_F(B4P_BallRingPairTest, calc_Sliding_1) {
	this->init_25BSWZ();

	int i = 0;	// 溝番号
	double a = 0.668610936573505e-3; // 接触楕円長径[m]
	double Pmax = 1;		// 最大面圧（使わない）
	double F_norm = 100;	// 転動体荷重[N]
	double fratio = 0.4;	// 油膜接触割合
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double sin_alp1 = 0.540517404112816;		// 外輪接触角の正弦(balac計算結果)
	double sin_alp2 = 0.540517404112816;		// 内輪接触角の正弦(balac計算結果)
	double cos_alp1 = 0.841332832980588;		// 外輪接触角の余弦(balac計算結果)
	double cos_alp2 = -0.841332832980588;		// 内輪接触角の余弦(balac計算結果)
	double dx_out = 4.389357066845456E-006;		// 外輪側接近量(balac計算結果)


	this->BOP.TR = new Tribology::Stab_Traction();	// トラクション係数を0.2に固定
	this->BOP.CL = new Tribology::Stab_Coulomb();	// クーロン摩擦係数を0.3に固定

	// 玉がY方向（周方向）に無回転で滑っている状態を仮定
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

	// 手計算の結果と比較
	double _fs = 100 * (0.4 * 0.2 + 0.6 * 0.3);
	double bp = this->BL.r - dx_out * 0.5;
	Vector3d _Fbs = Vector3d(0, -_fs, 0);
	Vector3d _Tbs = Vector3d(_fs *bp * cos_alp1, 0, -_fs * bp * sin_alp1);
	double err = 0.05;
	EXPECT_NEAR((Fbs - _Fbs).norm(), 0, _Fbs.norm()*err);
	EXPECT_NEAR((Tbs - _Tbs).norm(), 0, _Tbs.norm()*err);


	return;
}

// (2) 玉が外輪に対してスピン方向に回転した時の滑り摩擦を計算
// ただし，クーロン摩擦係数のみとした
TEST_F(B4P_BallRingPairTest, calc_Sliding_2) {
	this->init_25BSWZ();

	int i = 0;	// 溝番号
	double a = 0.668610936573505e-3; // 接触楕円長径[m]
	double Pmax = 1;		// 最大面圧（使わない）
	double F_norm = 100;	// 転動体荷重[N]
	double fratio = 0.0;	// 油膜接触割合
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double sin_alp1 = 0.540517404112816;		// 外輪接触角の正弦(balac計算結果)
	double cos_alp1 = 0.841332832980588;		// 外輪接触角の余弦(balac計算結果)
	double dx_out = 4.389357066845456E-006;		// 外輪側接近量(balac計算結果)


	this->BOP.TR = new Tribology::Stab_Traction();	// トラクション係数を0.2に固定
	this->BOP.CL = new Tribology::Stab_Coulomb();	// クーロン摩擦係数を0.3に固定

	// 玉が速度0で，接触面垂直方向を軸に回転している場合を仮定
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



	// 手計算の結果と比較
	double _fs = 100 * (0.4 * 0.2 + 0.6 * 0.3);
	double bp = this->BL.r - dx_out * 0.5;
	Vector3d _Fbs = Vector3d(0, 0, 0);
	EXPECT_NEAR((Fbs - _Fbs).norm(), 0, 0.1);			// 滑り摩擦のトータルが0
	Vector3d eg = Vector3d(sin_alp1, 0, cos_alp1);
	double dir = Tbs.dot(eg);
	Vector3d Dir = Tbs.cross(eg);
	EXPECT_LT(dir, 0);								// トルクの向きが接触面に対してボール側になっていることを確認
	EXPECT_NEAR(Dir.norm(), 0, 0.1);				// トルクが接触面に対して垂直であるか確認
	return;
}

// (3) 【外輪+側】玉接触角が30°でほぼ純転がりしているときの摩擦力の向きを評価
TEST_F(B4P_BallRingPairTest, calc_Sliding_3) {
	this->init_50BSWZ();
	int i = 1;							// 溝番号
	double a = 0.3160272198e-3;			// 接触楕円長径[m]
	double Pmax = 2e6;					// 最大面圧（使わない）
	double F_norm = 20;		// 転動体荷重[N]
	double fratio = 0.0;				// 油膜接触割合
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double sin_alp1 = 0.5;				// 30°のsin
	double cos_alp1 = sqrt(3) * 0.5;	// 30°のcos
	double dx_out = 0;		// 外輪側接近量(balac計算結果)

	// 流体潤滑：なし，クーロン摩擦：tanカーブ
	this->BOP.TR = new Tribology::Stab_Traction();
	this->BOP.CL = new Tribology::Tangent();
	this->BOP.FT = new Tribology::FilmThicknessNothing();
	// 玉が純転がり(時計回り)で，接触面垂直方向を軸に回転している場合を仮定
	double wx = 20;
	Vector3d bl_v = Vector3d(0, wx * (-dx_out + this->BL.r), 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d bl_w = Vector3d(wx * cos_alp1, 0, -wx * sin_alp1);

	// 玉位置は溝R中心の位置を基準に決定
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



	// 摩擦力・モーメントの評価
	Vector3d e_fs = Fbs.normalized();
	Vector3d e_ts = Tbs.normalized();
	Vector3d e_fsex(0, -1, 0);
	Vector3d e_tsex(cos_alp1, 0, -sin_alp1);
	EXPECT_NEAR((e_fs - e_fsex).norm(), 0, 1e-6);			// 滑り摩擦の向き
	EXPECT_NEAR((e_ts - e_tsex).norm(), 0, 1e-6);			// モーメントの向き

	return;
}


// (4) 【外輪+側】玉接触角が0°で玉自転軸がx軸基準で-45°を向いている場合の玉モーメントの向きが正しく計算できているか検証
TEST_F(B4P_BallRingPairTest, calc_Sliding_4) {
	this->init_50BSWZ();
	int i = 1;							// 溝番号
	double a = 0.3160272198e-3;			// 接触楕円長径[m]
	double Pmax = 2e6;					// 最大面圧（使わない）
	double F_norm = 20;					// 転動体荷重[N]
	double fratio = 0.0;				// 油膜接触割合
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double dx_out = 0;					// 外輪側接近量
	double cos_45 = 1.0 / sqrt(2.0); 
	double sin_45 = 1.0 / sqrt(2.0);

	this->BOP.TR = new Tribology::Stab_Traction();
	this->BOP.CL = new Tribology::Tangent();
	this->BOP.FT = new Tribology::FilmThicknessNothing();
	// 玉が速度0で，接触面垂直方向を軸に回転している場合を仮定
	double wx = 100;
	Vector3d bl_v = Vector3d(0, wx * (this->BL.r + dx_out * 0.5) * cos_45, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d bl_w = Vector3d(wx * sin_45, 0,  -wx * cos_45);

	// 玉位置は溝R中心の位置を基準に決定
	double R_out = this->OR.GV[i].r + dx_out - this->BL.r;
	Vector3d bl_x = Vector3d(this->OR.GV[i].Rx, 0, R_out + this->OR.GV[i].Rr);
	this->BL.set_y(bl_x, bl_v, q, bl_w);
	Vector3d zero = Vector3d::Zero();
	this->OR.set_y(zero, zero, q, zero);
	double rr = this->OR.GV[i].r + dx_out * 0.5;
	Vector3d p = Vector3d(this->OR.GV[i].Rx, 0, rr + this->OR.GV[i].Rr);


	this->BOP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);

	// 手計算の結果と比較
	EXPECT_GT(Tbs.z(), 0);								// z+方向にモーメントが働く
	EXPECT_LT(Tis.z(), 0);								// z-方向にモーメントが働く

	return;
}

// (5) 【外輪+側】玉接触角が0°で玉自転軸がx軸基準で+45°を向いている場合の玉モーメントの向きが正しく計算できているか検証
TEST_F(B4P_BallRingPairTest, calc_Sliding_5) {
	this->init_50BSWZ();
	int i = 1;							// 溝番号
	double a = 0.3160272198e-3;			// 接触楕円長径[m]
	double Pmax = 2e6;					// 最大面圧（使わない）
	double F_norm = 20;		// 転動体荷重[N]
	double fratio = 0.0;				// 油膜接触割合
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double dx_out = 0;		// 外輪側接近量
	double cos_45 = 1.0 / sqrt(2.0);
	double sin_45 = 1.0 / sqrt(2.0);

	this->BOP.TR = new Tribology::Stab_Traction();
	this->BOP.CL = new Tribology::Tangent();
	this->BOP.FT = new Tribology::FilmThicknessNothing();
	// 玉が速度0で，接触面垂直方向を軸に回転している場合を仮定
	double wx = 100;
	Vector3d bl_v = Vector3d(0, wx * (this->BL.r + dx_out * 0.5) * cos_45, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d bl_w = Vector3d(wx * sin_45, 0, wx * cos_45);

	// 玉位置は溝R中心の位置を基準に決定
	double R_out = this->OR.GV[i].r + dx_out - this->BL.r;
	Vector3d bl_x = Vector3d(this->OR.GV[i].Rx, 0, R_out + this->OR.GV[i].Rr);
	this->BL.set_y(bl_x, bl_v, q, bl_w);
	Vector3d zero = Vector3d::Zero();
	this->OR.set_y(zero, zero, q, zero);
	double rr = this->OR.GV[i].r + dx_out * 0.5;
	Vector3d p = Vector3d(this->OR.GV[i].Rx, 0, rr + this->OR.GV[i].Rr);


	this->BOP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);
	// 手計算の結果と比較
	EXPECT_LT(Tbs.z(), 0);								// z-方向にモーメントが働く

	return;
}

// (6) 【内輪-側】玉接触角が0°で玉自転軸がx軸基準で+45°を向いている場合の玉モーメントの向きが正しく計算できているか検証
TEST_F(B4P_BallRingPairTest, calc_Sliding_6) {
	this->init_50BSWZ();
	int i = 0;							// 溝番号
	double a = 0.3160272198e-3;			// 接触楕円長径[m]
	double Pmax = 2e6;					// 最大面圧（使わない）
	double F_norm = 20;		// 転動体荷重[N]
	double fratio = 0.0;				// 油膜接触割合
	Vector3d Fbs;
	Vector3d Tbs, Tis;
	double dxi = 0;		// 内輪側接近量
	double cos_45 = 1.0 / sqrt(2.0);
	double sin_45 = 1.0 / sqrt(2.0);
	double cos_m135 = -1.0 / sqrt(2.0);
	double sin_m135 = -1.0 / sqrt(2.0);
	this->BIP.TR = new Tribology::Stab_Traction();
	this->BIP.CL = new Tribology::Tangent();
	this->BIP.FT = new Tribology::FilmThicknessNothing();
	// 接触点における滑り速度が0の時を仮定
	double wx = 100;
	Vector3d bl_v = Vector3d(0, wx * (this->BL.r + dxi * 0.5) * cos_45, 0);
	Quaterniond q = Quaterniond(1, 0, 0, 0);
	Vector3d bl_w = Vector3d(wx * sin_45, 0, wx * cos_45);

	// 玉位置は溝R中心の位置を基準に決定
	double R_in = this->IR.GV[i].r + dxi - this->BL.r;
	Vector3d bl_x = Vector3d(this->IR.GV[i].Rx, 0, this->IR.GV[i].Rr - R_in);
	this->BL.set_y(bl_x, bl_v, q, bl_w);
	Vector3d zero = Vector3d::Zero();
	this->IR.set_y(zero, zero, q, zero);
	double rr = this->IR.GV[i].r + dxi * 0.5;
	Vector3d p = Vector3d(this->IR.GV[i].Rx, 0, rr + this->IR.GV[i].Rr);
	this->BIP.calc_Sliding(i, p, a, Pmax, F_norm, fratio, Fbs, Tbs, Tis);
	EXPECT_LT(Tbs.z(), 0);								// z-方向にモーメントが働く

	return;
}


// (1) 正の接触角の時の外輪接触角を計算し．手計算の数値と比較
TEST_F(B4P_BallRingPairTest, ContactAngle_1) {
	this->init_25BSWZ();

	double err = 0.01;
	double sin_alp1 = 0.5;				// 外輪接触角の正弦
	double cos_alp1 = 0.8660254;		// 外輪接触角の余弦
	double alp = this->BOP.ContactAngle(cos_alp1, sin_alp1);
	double _alp = Unit::deg2rad(30);
	EXPECT_NEAR(_alp, alp, alp * err);
	return;
}

// (2) 負の接触角の時の外輪接触角を計算し．手計算の数値と比較
TEST_F(B4P_BallRingPairTest, ContactAngle_2) {
	this->init_25BSWZ();

	double err = 0.01;
	double sin_alp1 = -0.5;				// 外輪接触角の正弦
	double cos_alp1 = 0.8660254;		// 外輪接触角の余弦
	double alp = this->BOP.ContactAngle(cos_alp1, sin_alp1);
	double _alp = Unit::deg2rad(-30);
	EXPECT_NEAR(_alp, alp, abs(alp) * err);
	return;
}

// (3) 正の接触角の時の内輪接触角を計算し．手計算の数値と比較
TEST_F(B4P_BallRingPairTest, ContactAngle_3) {
	this->init_25BSWZ();

	double err = 0.01;
	double sin_alp1 = 0.5;				// 外輪接触角の正弦
	double cos_alp1 = -0.8660254;		// 外輪接触角の余弦
	double alp = this->BIP.ContactAngle(cos_alp1, sin_alp1);
	double _alp = Unit::deg2rad(30);
	EXPECT_NEAR(_alp, alp, alp * err);
	return;
}

// (4) 負の接触角の時の内輪接触角を計算し．手計算の数値と比較
TEST_F(B4P_BallRingPairTest, ContactAngle_4) {
	this->init_25BSWZ();

	double err = 0.01;
	double sin_alp1 = -0.5;				// 外輪接触角の正弦
	double cos_alp1 = -0.8660254;		// 外輪接触角の余弦
	double alp = this->BIP.ContactAngle(cos_alp1, sin_alp1);
	double _alp = Unit::deg2rad(-30);
	EXPECT_NEAR(_alp, alp, abs(alp) * err);
	double alp__ = Unit::rad2deg(atan2(sin_alp1, cos_alp1));
	return;
}

// (1) 条件分岐の確認
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

