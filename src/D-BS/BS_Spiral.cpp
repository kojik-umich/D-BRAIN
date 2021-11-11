/*******************************************************************************
!								"BS_Spiral.cpp"
!													2020/01/07	[Core-T]	楠崎
!	Siralを継承して作成しているクラス．
!	ゴシックアーク溝をもつため，追加のメンバ変数としてGroove構造体をもつ．
!
!*******************************************************************************/
#include "BS_Spiral.h"


void BS_Spiral::Nothing(void) {
}

void BS_Spiral::init(const BS_In::Cylinder::Spiral&spiral) {

	this->SP.init(spiral.alp, spiral.l, spiral.r);
	double nd = this->SP.get_nd();
	double r = this->SP.get_r();

	for (int ig = 0; ig < 2; ig++) {
		this->GV[ig].eta   = Vector2d(spiral.groove[ig].eta);
		this->GV[ig].r     = spiral.groove[ig].r;
		this->GV[ig].r_inv = 1.0 / spiral.groove[ig].r;
		this->GV[ig].sigma = spiral.groove[ig].sigma;
		this->GV[ig].R     = nd * nd / r - this->GV[ig].eta[0];	// ネジの捩率を考慮した R = nd^2 / r
	}
	return;
}

// ボールの接触角の余弦から，螺旋の曲率2成分を返すメソッド．
Vector2d BS_Spiral::get_rho
(					// out:	[1/m]	曲率．第0成分：螺旋方向．第1成分：溝曲率方向．
	double cos_alp,	// in:	[-]		接触角のcos．（ナットラジアル当たり=0度，シャフトラジアル当たり=180度）
	int ig			// in:	[-]		溝番号．
) {
	double R = this->GV[ig].R;

	Vector2d rho = Vector2d(
		cos_alp / (R - this->GV[ig].r * cos_alp),
		-this->GV[ig].r_inv
	);
	return rho;
}

