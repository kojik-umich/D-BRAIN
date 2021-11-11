/*******************************************************************************
!								"BS_Spiral.cpp"
!													2020/01/07	[Core-T]	���
!	Siral���p�����č쐬���Ă���N���X�D
!	�S�V�b�N�A�[�N�a�������߁C�ǉ��̃����o�ϐ��Ƃ���Groove�\���̂����D
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
		this->GV[ig].R     = nd * nd / r - this->GV[ig].eta[0];	// �l�W�̝������l������ R = nd^2 / r
	}
	return;
}

// �{�[���̐ڐG�p�̗]������C�����̋ȗ�2������Ԃ����\�b�h�D
Vector2d BS_Spiral::get_rho
(					// out:	[1/m]	�ȗ��D��0�����F���������D��1�����F�a�ȗ������D
	double cos_alp,	// in:	[-]		�ڐG�p��cos�D�i�i�b�g���W�A��������=0�x�C�V���t�g���W�A��������=180�x�j
	int ig			// in:	[-]		�a�ԍ��D
) {
	double R = this->GV[ig].R;

	Vector2d rho = Vector2d(
		cos_alp / (R - this->GV[ig].r * cos_alp),
		-this->GV[ig].r_inv
	);
	return rho;
}

