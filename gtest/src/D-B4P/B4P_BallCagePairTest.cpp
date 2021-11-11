#include "pch.h"
#include "bal_StabCage.h"
#include "B4P_StabIn.h"
#include "B4P_StabOut.h"

/*************************************************************************************************
BallCageRingPairクラスのテスト基本方針
下位クラスで玉の接触位置や接近量については検証できているため，このクラスでは上記については検証せず，
接近量や接触位置が分かっているときの荷重が正しく計算できていることだけを検証する．
また，保持器クラスはスタブを用いる．
**************************************************************************************************/

class B4P_BallCagePairTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	B4P_BallCagePair  BCP;
	Ball*BL;
	bal_Cage*CG;
	B4P_StabOut OUT;
	Vector3d BL_Cent;	// ポケット中心

	virtual void SetUp() {
		B4P_StabIn FI;
		this->BL = new Ball;
		int msmax = 21;
		OUT.allocate(1, _MAX_CONTACT_, msmax);
		this->CG = new bal_StabCage;
		this->CG->x = Vector3d::Zero();
		this->CG->v = Vector3d::Zero();
		this->CG->w = Vector3d::Zero();
		this->CG->q = Quaterniond::Identity();
		this->CG->m = 0.2;
		this->BL_Cent = Vector3d(0.0, 0.0, 30.0e-3);
		double D = 10.0e-3;				// 玉径[m]
		double E = 207760000000;		// ヤング率[Pa]
		double por = 0.29;				// ポアソン比[-]
		double den = 7830;				// 密度[kg/m^3]
		double rms = 0.00000002;		// 粗さrms[m]
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}
		this->BL->init(D, E, por, den, rms, x_const, Rx_const);
		this->BL->x = Vector3d::Zero();
		this->BL->v = Vector3d::Zero();
		this->BL->w = Vector3d::Zero();
		this->BL->q = Quaterniond::Identity();
		this->BL->m = 0.2;
		this->BCP.link(this->BL, this->CG, 0);
		FI.BCP.mu = 0.3;
		FI.BCP.dzeta = 0.2;
		this->BCP.init(FI);
	}

	virtual void TearDown() {
	}
};

// 接近量，接近位置が分かっているときの荷重を計算
// calc_forceとmake_OutParamをセットでテスト
TEST_F(B4P_BallCagePairTest, calc_force_test1) {
	// 本来ならば get_ContactPoint で接触点候補・剛性を計算するが，
	// get_ContactPoint で行っている計算は検証対象としないため，
	// スタブを用いて接触位置が(0, 5e-3, 30.0e-3)となるように固定
	Vector3d dx = Vector3d(0, 1e-3, 0);	// 玉と保持器ポケットの接近量
	this->BL->x = this->BL_Cent + dx;
	this->BL->v = Vector3d(0, 1, 0);
	this->BL->w = Vector3d(1, 0, 0);
	Vector3d Fbc, Tbc, Fcb, Tcb; // Fbc : ボール(b)が保持器(c)から受ける力．
	this->BCP.calc_force(Fbc, Tbc, Fcb, Tcb);

	// スタブ保持器を用いているため，接触位置は固定 → 接近量は 1e-3 m，接近速度は 1 m/s
	// 接触剛性は 1e5 N/m, 減衰係数は 40 → 弾性接触荷重は (0, -140, 0)
	// 摩擦力は(0, 0, -42)
	Vector3d Fbc_ex = Vector3d(0, -140, -42);
	Vector3d Fcb_ex = Vector3d(0, 140, 42);
	Vector3d Tbc_ex = Vector3d(-0.168, 0, 0);
	Vector3d Tcb_ex = Vector3d(-3.99, 0, 0);
	EXPECT_LT((Fbc_ex - Fbc).norm(), 0.1);
	EXPECT_LT((Fcb_ex - Fcb).norm(), 0.1);
	EXPECT_LT((Tbc_ex - Tbc).norm(), 0.1);
	EXPECT_LT((Tcb_ex - Tcb).norm(), 0.1);

	this->BCP.save(OUT.BCP[0]);
	Vector3d p_mm(OUT.BCP[0].CP[1].p);
	Vector3d p_mm_ex = Vector3d(0, 5.0e-3, 30.0e-3);
	EXPECT_LT((p_mm_ex - p_mm).norm(), 0.1);
	Vector3d Fn_(OUT.BCP[0].CP[1].Fn_);
	Vector3d Fs_(OUT.BCP[0].CP[1].Fs_);
	Vector3d F_ex = Vector3d(0, -140, -42);
	EXPECT_LT((F_ex - Fn_ - Fs_).norm(), 0.1);

	return;
};

// (1)すべり速度のテスト
TEST_F(B4P_BallCagePairTest, get_us_1) {
	Vector3d dx = Vector3d(0, 1e-3, 0);	// 玉と保持器ポケットの接近量
	this->BL->x = this->BL_Cent + dx;
	this->BL->v = Vector3d(0, 1, 0);
	this->BL->w = Vector3d(1, 0, 0);
	Vector3d p(0, 5e-3, 30.0e-3);
	Vector3d us = this->BCP.get_us(p);
	Vector3d us_ans(0, 0, 4e-3);
	EXPECT_NEAR((us - us_ans).norm(), 0, 1e-3);
	return;
}