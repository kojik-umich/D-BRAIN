#define EXPORT

#include "ball.h"
#include "BS_Nut.h"
#include "BS_Shaft.h"
#include "BS_BallCylinderPair.h"
#include "BS_SimulinkIn.h"
#include "BS_Simulink.h"

BS_Simulink *newSimulink() {

	BS_Simulink*Simulink = new BS_Simulink();

	Simulink->BS_BallNutPair0   = new BS_BallCylinderPair();
	Simulink->BS_BallNutPair1   = new BS_BallCylinderPair();
	Simulink->BS_BallShaftPair0 = new BS_BallCylinderPair();
	Simulink->BS_BallShaftPair1 = new BS_BallCylinderPair();
	Simulink->BS_SimulinkIn     = new BS_SimulinkIn();

	return Simulink;
}

void BS_read(BS_Simulink*Simulink, const char dbsin_csv[]) {

	BS_SimulinkIn *SLIN = static_cast<BS_SimulinkIn *>(Simulink->BS_SimulinkIn);

	SLIN->read_csv(dbsin_csv);

	return;
}

void BS_init(BS_Simulink*Simulink) {

	BS_SimulinkIn *SLIN = static_cast<BS_SimulinkIn *>(Simulink->BS_SimulinkIn);
	BS_BallCylinderPair *BNP0 = static_cast<BS_BallCylinderPair *>(Simulink->BS_BallNutPair0);
	BS_BallCylinderPair *BNP1 = static_cast<BS_BallCylinderPair *>(Simulink->BS_BallNutPair1);

	BNP0->init(SLIN->BallNutPair[0], SLIN->tribology, SLIN->oil);
	BNP1->init(SLIN->BallNutPair[1], SLIN->tribology, SLIN->oil);

	return;
}

void deleteBS(BS_Simulink*Simulink) {
	delete Simulink;
}




