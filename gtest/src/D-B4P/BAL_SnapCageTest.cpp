#include "pch.h"
#include "B4P_StabIn.h"

#define _FACE_		1e5
#define _CORNEROUT_ 2e5
#define _OPEN_		3e5
#define _CORNERIN_	4e5
#define _EDGEIN_	5e5
#define _EDGEOUT_	6e5

class bal_SnapCageTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	bal_SnapCage CG;
	double BL_r;
	Vector3d BL_Cent;						

	virtual void SetUp() {
		// 物性値は50BSWZのものを用いる．
		B4P_StabIn FI;

		FI.Snap.dout		= 72.2e-3;
		FI.Snap.din			= 62.8e-3;
		FI.Cage.m			= 65.553e-3;

		Vector3d rmg0(3.64545273e-3, 0.0, 0.0);
		for (int i = 0; i < 3; i++)
			FI.Cage.rmg0[i]		= rmg0[i];

		FI.ballnum			= 11;
		FI.Cage.Ix			= 7.5031;
		FI.Cage.Iyz			= 0.022;

		FI.balldia			= 11.1e-3;		// 玉径[m]
		FI.Snap.R			= 5.74e-3;		// ポケット半径[m]
		FI.Snap.ropen		= 5.055e-3;		// 開口部半径[m]
		FI.Snap.Kface = _FACE_;
		FI.Snap.Kcornerout = _CORNEROUT_;
		FI.Snap.Kopen = _OPEN_;
		FI.Snap.Kcornerin = _CORNERIN_;
		FI.Snap.Kedgein = _EDGEIN_;
		FI.Snap.Kedgeout = _EDGEOUT_;
		FI.Snap.jc			= 0.0;

		this->CG.init(FI);

		// 保持器幾何中心が原点に来るように調整．
		for (int i = 0; i < 3; i++)
			this->CG.x[i] = -FI.Cage.rmg0[i];
		this->CG.q = Quaterniond::Identity();

		this->BL_r = FI.balldia / 2;
		this->BL_Cent = Vector3d(0.0, 0.0, 67.5e-3 / 2);
	}

	virtual void TearDown() {
	}
};

// (1) 角座標，垂直高さを検証
TEST_F(bal_SnapCageTest, init0) {
	// ポケット0の座標
	Vector3d corno0_ex((3.5 - 0.78070335) * 1e-3, -4.61921891e-3, 35.80325148e-3);
	Vector3d corno1_ex((3.5 - 0.78070335) * 1e-3, 4.61921891e-3, 35.80325148e-3);
	Vector3d corni0_ex((3.5 - 0.78070335) * 1e-3, -4.30670957e-3, 31.10325148e-3);
	Vector3d corni1_ex((3.5 - 0.78070335) * 1e-3, 4.30670957e-3, 31.10325148e-3);

	// ポケット頂点の座標は回転行列を用いて計算し，プログラム計算値と比較
	VectorXd th = VectorXd(12).setLinSpaced(0, 2 * Numeric::pi);
	for (int i = 0; i < 11; i++) {
		double thi = th[i];
		Matrix3d roti;		// 回転行列
		roti << 1, 0, 0,
			0, cos(-th[i]), sin(-th[i]),
			0, -sin(-th[i]), cos(-th[i]);
		EXPECT_NEAR((this->CG.PK[i].corno0 - roti * corno0_ex).norm(), 0, 1e-6);
		EXPECT_NEAR((this->CG.PK[i].corno1 - roti * corno1_ex).norm(), 0, 1e-6);
		EXPECT_NEAR((this->CG.PK[i].corni0 - roti * corni0_ex).norm(), 0, 1e-6);
		EXPECT_NEAR((this->CG.PK[i].corni1 - roti * corni1_ex).norm(), 0, 1e-6);
	}


	// 最小垂直高さと開き角α(開口部高さで評価), γはすべてのポケットで同じ数値になっていることを確認．
	for (int i = 0; i < 11; i++){
		EXPECT_NEAR(this->CG.PK[i].zomin, 35.69370222e-3, 1e-6);
		EXPECT_NEAR(this->CG.PK[i].zimin, 30.99370222e-3, 1e-6);
		double h0_ex = 3.5e-3 - 0.78070681e-3;
		EXPECT_NEAR(this->CG.PK[i].h0, h0_ex, 1e-6);

		// cosγは垂直高さを用いて定義通り計算
		double cos_gammao_ex = (35.80325148e-3 - 33.75e-3) / 5.055e-3;
		double cos_gammai_ex = (31.10325148e-3 - 33.75e-3) / 5.055e-3;
		
		EXPECT_NEAR(this->CG.PK[i].cos_gammao, cos_gammao_ex, 1e-6);
		EXPECT_NEAR(this->CG.PK[i].cos_gammai, cos_gammai_ex, 1e-6);
	}
	return;
}

// (1) 玉座標を入力，出力される接触点の数と向きを評価
// 関数 get_ContactPoint・how_Contact・where_Contact がテスト対象
TEST_F(bal_SnapCageTest, get_ContactPoint1) {

	Vector3d BL_x;
	Vector3d return_[_MAX_CONTACT_];
	double k[_MAX_CONTACT_];
	Vector3d dx;
	int nc, ptt[_MAX_CONTACT_];
	Vector3d _pb, e_pb, _pp, e_pp;

	// ケース(A) face 
	dx = Vector3d(0.0, 0.4e-3, 0.0);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	double dot = dx.dot(return_[0] - BL_x) / dx.norm() / (return_[0] - BL_x).norm();
	EXPECT_EQ(nc, 1);										// 接触点の数が1つ
	EXPECT_NEAR(dot, 1, 1e-6);								// 玉のポケット内方向ベクトルと接近量ベクトルが同じ向き
	double La = (return_[0] - this->CG.PK[0].x).norm();
	EXPECT_NEAR(La, this->CG.PK[0].R, 1e-6);				// ポケット球上に接触点が存在
	EXPECT_NEAR(k[0], _FACE_, 1e-3);							// 接触剛性

	// ケース(B) 外輪エッジ当たり edgeo
	dx = Vector3d(0.0, 0.4e-3, 0.4e-3);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_EQ(nc, 1);										// 接触点の数が1つ
	double Lb1 = return_[0].norm();
	EXPECT_NEAR(Lb1, this->CG.ro, 1e-6);					// 保持器外周上に点が乗っているか確認
	double Lb2 = (return_[0] - this->CG.PK[0].x).norm();
	EXPECT_NEAR(Lb2, this->CG.PK[0].R, 1e-6);				// ポケット球に点が乗っているか確認

	// 径方向(Z軸方向)からみたときの方位ベクトルが玉の方位ベクトルと同じであるか確認
	_pb = return_[0] - this->CG.PK[0].x;
	e_pb = Vector3d(_pb[0], _pb[1], 0).normalized();		// 径方向からみた接触点の方位ベクトル(ポケット中心基準)				
	_pp = return_[0] - this->CG.PK[0].x;
	e_pp = Vector3d(_pp[0], _pp[1], 0).normalized();		// 径方向からみた玉の方位ベクトル(ポケット中心基準)
	EXPECT_NEAR((e_pb - e_pp).norm(), 0, 1e-3);
	EXPECT_NEAR(k[0], _EDGEOUT_, 1e-3);							// 接触剛性
		
	// ケース(C) 内輪エッジ当たり edgei
	dx = Vector3d(0.0, 0.4e-3, -0.4e-3);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_EQ(nc, 1);										// 接触点の数が1つ
	double Lc1 = return_[0].norm();
	EXPECT_NEAR(Lc1, this->CG.ri, 1e-6);					// 保持器内側の円周上に点が乗っているか確認
	double Lc2 = (return_[0] - this->CG.PK[0].x).norm();
	EXPECT_NEAR(Lc2, this->CG.PK[0].R, 1e-6);				// ポケット球に点が乗っているか確認

	// 径方向(Z軸方向)からみたときの方位ベクトルが玉の方位ベクトルと同じであるか確認
	_pb = return_[0] - this->CG.PK[0].x;
	e_pb = Vector3d(_pb[0], _pb[1], 0).normalized();		// 径方向からみた接触点の方位ベクトル(ポケット中心基準)				
	_pp = return_[0] - this->CG.PK[0].x;
	e_pp = Vector3d(_pp[0], _pp[1], 0).normalized();		// 径方向からみた玉の方位ベクトル(ポケット中心基準)
	EXPECT_NEAR((e_pb - e_pp).norm(), 0, 1e-3);
	EXPECT_NEAR(k[0], _EDGEIN_, 1e-3);							// 接触剛性

	// ケース(D) 開口部当たり aperture
	dx = Vector3d(0.4e-3, 0.4e-3, 0.0);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_EQ(nc, 2);										// 接触点の数が2つ
	EXPECT_NEAR(return_[0].x(), this->CG.PK[0].h0, 1);		// 上面に点が乗っているか確認
	double Ld1 = (return_[0] - this->CG.PK[0].x).norm();
	EXPECT_NEAR(Ld1, this->CG.PK[0].R, 1e-6);				// ポケット球に点が乗っているか確認

	// 上面(X軸方向)からみたときの接触点の方位ベクトルが玉の方位ベクトルと同じであるか確認（方向ベクトルはポケット中心基準に計算）
	_pb = return_[0] - this->CG.PK[0].x;
	e_pb = Vector3d(0, _pb[1], _pb[2]).normalized();		// 上面からみた接触点の方位ベクトル(ポケット中心基準)				
	_pp = return_[0] - this->CG.PK[0].x;
	e_pp = Vector3d(0, _pp[1], _pp[2]).normalized();		// 上面からみた玉の方位ベクトル(ポケット中心基準)
	EXPECT_NEAR((e_pb - e_pp).norm(), 0, 1e-3);
	EXPECT_NEAR(k[0], _OPEN_, 1e-3);							// 接触剛性
	
	// ケース(E) 角外輪側 cornero
	dx = Vector3d(0.4e-3, -0.4e-3, 0.4e-3);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR((return_[0] - this->CG.PK[0].corno0).norm(), 0, 1e-3);
	EXPECT_NEAR((return_[1] - this->CG.PK[0].corno1).norm(), 0, 1e-3);		// 接触点が頂点と一致
	EXPECT_EQ(nc, 2);														// 接触点の数が2つ
	EXPECT_NEAR(k[0], _CORNEROUT_, 1e-3);											// 接触剛性

	// ケース(F) 角内輪側 corneri
	dx = Vector3d(0.4e-3, -0.4e-3, -0.4e-3);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR((return_[0] - this->CG.PK[0].corni0).norm(), 0, 1e-3);
	EXPECT_NEAR((return_[1] - this->CG.PK[0].corni1).norm(), 0, 1e-3);		// 接触点が頂点と一致
	EXPECT_EQ(nc, 2);												  		// 接触点の数が2つ
	EXPECT_NEAR(k[0], _CORNERIN_, 1e-3);											// 接触剛性

	// ケース(O) exception（初期位置）
	dx = Vector3d(0.0, 0.0, 0.0);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_EQ(nc, 0);			// 接触点なし

	// ケース(I) exception（保持器幾何中心・ポケット中心・ボール中心が一列）
	dx = Vector3d(0.0, 0.0, 1e-3);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_EQ(nc, 0);			// 接触点なし

	// ケース(J) exception（軸とポケット中心-ボール中心が並行）
	dx = Vector3d(1e-3, 0.0, 0.0);
	BL_x = this->BL_Cent + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_EQ(nc, 0);			// 接触点なし

	return;
};

// (2) 頂点の境界値付近の当たり判定の確認（接触点の座標はテストしない）
TEST_F(bal_SnapCageTest, get_ContactPoint2) {

	Vector3d BL_x;
	Vector3d return_[_MAX_CONTACT_];
	double k[_MAX_CONTACT_];
	Vector3d dx;
	double dot;
	int nc, ptt[_MAX_CONTACT_];

	// 角の位置に玉位置を設定(外輪側)
	Vector3d corno0_ex((3.5 - 0.78070335) * 1e-3, -4.61921891e-3, 35.80325148e-3);

	// -x, +zであればエッジ当たり
	dx = Vector3d(-1e-6, 0, 1e-6);
	BL_x = corno0_ex + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _EDGEOUT_, 1e-3);			// 接触剛性
	EXPECT_EQ(nc, 1);						// 接触点の数が1つ

	// +x, +zであれば角当たり
	dx = Vector3d(1e-6, 0, 1e-6);
	BL_x = corno0_ex + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _CORNEROUT_, 1e-3);			// 接触剛性
	EXPECT_EQ(nc, 2);						// 接触点の数が2つ

	// +x, -zであれば開口部当たり
	dx = Vector3d(1e-6, 0, -1e-6);
	BL_x = corno0_ex + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _OPEN_, 1e-3);			// 接触剛性
	EXPECT_EQ(nc, 2);						// 接触点の数が2つ


	// 角の位置に玉位置を設定(内輪側)
	Vector3d corni0_ex((3.5 - 0.78070335) * 1e-3, -4.30670957e-3, 31.10325148e-3);

	// -x, -zであればエッジ当たり
	dx = Vector3d(-1e-6, 0, -1e-6);
	BL_x = corni0_ex + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _EDGEIN_, 1e-3);			// 接触剛性
	EXPECT_EQ(nc, 1);						// 接触点の数が1つ

	// +x, -zであれば角当たり
	dx = Vector3d(1e-6, 0, -1e-6);
	BL_x = corni0_ex + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _CORNERIN_, 1e-3);			// 接触剛性
	EXPECT_EQ(nc, 2);						// 接触点の数が2つ

	// +x, +zであれば開口部当たり
	dx = Vector3d(1e-6, 0, 1e-6);
	BL_x = corni0_ex + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _OPEN_, 1e-3);			// 接触剛性
	EXPECT_EQ(nc, 2);						// 接触点の数が2つ

	return;
}

// (3) 開口部の境界値付近の当たり判定の確認（接触点の座標はテストしない）
TEST_F(bal_SnapCageTest, get_ContactPoint3) {
	
	Vector3d BL_x;
	Vector3d return_[_MAX_CONTACT_];
	double k[_MAX_CONTACT_];
	Vector3d dx;
	double dot;
	int nc, ptt[_MAX_CONTACT_];

	// 開口部の座標
	Vector3d _x = Vector3d(this->CG.PK[0].h0, this->CG.PK[0].ropen, 0) 
		+ this->CG.PK[0].x;

	// -x方向に動くと面当たり
	dx = Vector3d(-1e-6, 0, 0);
	BL_x = _x + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _FACE_, 1e-3);			// 接触剛性
	EXPECT_EQ(nc, 1);						// 接触点の数が1つ

	// +x方向に動くと開口部当たり
	dx = Vector3d(1e-6, 0, 0);
	BL_x = _x + dx;
	nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, return_, k, ptt);
	EXPECT_NEAR(k[0], _OPEN_, 1e-3);			// 接触剛性
	EXPECT_EQ(nc, 2);						// 接触点の数が2つ
	
	return;
}

// (4) 境界値を跨ぐ付近（面→エッジ→コーナー→開口部→面）の接触点の移動を見，この近傍で連続になるかの確認．
TEST_F(bal_SnapCageTest, get_ContactPoint4) {

	double r = 1e-4;
	int n = 50;
	VectorXd th = VectorXd(n).setLinSpaced(0, 2* Numeric::pi);
	double dth = abs(th[1] - th[0]);
	double dx = r * dth;	// およそ 1.282e-05 m ⇒ 12 um

	// まずは，内コーナー周りで円形に玉位置を動かし，接触点の移動を確認する．．
	Vector3d corni0 = this->CG.PK[0].corni0;
	Vector3d ci0 = this->CG.PK[0].corni0 - this->CG.PK[0].x;
	Vector3d axi = Vector3d::UnitX().cross(ci0).normalized();
	Vector3d ayi = axi.cross(ci0).normalized();
	int ptt[_MAX_CONTACT_];
	Vector3d xc0;
	for (size_t i = 0; i < n; i++) {
		Vector3d BL_x = corni0 + r * cos(th[i]) * axi + r * sin(th[i]) * ayi;

		Vector3d xc[_MAX_CONTACT_];		double k[_MAX_CONTACT_];
		int nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, xc, k, ptt);
		if (i > 0)
			EXPECT_LT((xc0 - xc[0]).norm(), 3 * dx);	// ここを1*dxに変えるといくつかテストを通らないものが出てきます．下のコメントとセットでやるとより分かりやすいです．試して，ちょっと考えてみてみてください笑
		xc0 = xc[0];
		// cout << k[0] << endl;	// 剛性を出力させることで，接触判定がどうなっているのか確認できます．動作確認時にはこのコメントを外して実行してください．
	}

	// 同様に，外コーナー周り．
	Vector3d corno0 = this->CG.PK[0].corno0;
	Vector3d co0 = this->CG.PK[0].corno0 - this->CG.PK[0].x;
	Vector3d axo = Vector3d::UnitX().cross(co0).normalized();
	Vector3d ayo = axo.cross(co0).normalized();

	for (size_t i = 0; i < n; i++) {
		Vector3d BL_x = corno0 + r * cos(th[i]) * axo + r * sin(th[i]) * ayo;

		Vector3d xc[_MAX_CONTACT_];		double k[_MAX_CONTACT_];
		int nc = this->CG.get_ContactPoint(BL_x, this->BL_r, 0, xc, k, ptt);
		if (i > 0)
			EXPECT_LT((xc0 - xc[0]).norm(), 3 * dx);
		xc0 = xc[0];
		// cout << k[0] << endl;	// 剛性を出力させることで，接触判定がどうなっているのか確認できます．動作確認時にはこのコメントを外して実行してください．
	}
	return;
}

