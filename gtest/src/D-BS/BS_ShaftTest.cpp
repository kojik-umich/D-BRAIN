#include "pch.h"
#include "BS_FileIn_stab.h"

class BS_ShaftTest : public ::testing::Test {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

protected:
	BS_FileIn_stab IN;
	BS_Shaft ST;
	
	virtual void SetUp() {
	}
	void init(){
		this->IN.shaft.resize(1);
		this->IN.shaft[0].density = 1.0;
		this->IN.shaft[0].poisson = 1.0;
		this->IN.shaft[0].young = 1.0;
		this->IN.shaft[0].ri = 1.2;
		this->IN.shaft[0].ro = 1.2;
		
		this->IN.shaft[0].spiral.resize(1);

		this->IN.shaft[0].spiral[0].alp = 0.0;
		this->IN.shaft[0].spiral[0].l = 0.8;
		this->IN.shaft[0].spiral[0].r   = 1.0;
		
		for (int i = 0; i < 2; i++) {
			this->IN.shaft[0].spiral[0].groove[i].eta[0] = 0.0;
			this->IN.shaft[0].spiral[0].groove[i].eta[1] = 0.0;
			this->IN.shaft[0].spiral[0].groove[i].r = 1.0e-3;
			this->IN.shaft[0].spiral[0].groove[i].sigma = 1.0e-6;
		}

		this->ST.allocate(this->IN.shaft);

		double v0 = 0.0;
		double w0 = 0.0;
		bool x_const[3], Rx_const[3];
		for (int i = 0; i < 3; i++) {
			x_const[i] = false;
			Rx_const[i] = false;
		}
		this->ST.init(this->IN.shaft, x_const, Rx_const, v0, w0);

		Rigid::l = 1e-2;
		Rigid::t = 1e-3;

		this->ST.init_pos(v0, w0);

		this->ST.x = Vector3d(0.0, 1.0, 0.0);
		this->ST.set_dx();
	}

	virtual void TearDown() {
	}
};

TEST_F(BS_ShaftTest, to_etacoord0) {
	this->init();
	double th = 1e1;

	Vector3d eta0(th, 1e-3, 2e-3);
	Vector3d x0 = this->ST.CY->to_inertialcoord(0, eta0);

	Vector3d eta1 = this->ST.CY->to_etacoord(0, x0);

	double delta = (eta0 - eta1).norm();

	EXPECT_GT(1e-2, delta);
};
