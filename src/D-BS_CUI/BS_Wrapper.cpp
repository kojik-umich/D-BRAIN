#define EXPORT

#include "BS_BallScrew.h"
#include "BS_FileIn.h"
#include "BS_FileOut.h"
#include "BS_Wrapper.h"

DynamicBallscrew *newBS() {

	DynamicBallscrew*DBS = new DynamicBallscrew();

	DBS->BS_BallScrew = new BS_BallScrew();
	DBS->BS_FileIn    = new BS_FileIn();
	DBS->BS_FileOut   = new BS_FileOut();

	return DBS;
}

void BS_read(DynamicBallscrew*DBS, const char dbsin_csv[]) {

	BS_FileIn *FI = static_cast<BS_FileIn *>(DBS->BS_FileIn);

	FI->read_input_all(dbsin_csv);

	return;
}

void BS_init(DynamicBallscrew*DBS) {

	BS_FileIn *FI = static_cast<BS_FileIn *>(DBS->BS_FileIn);
	BS_FileOut *FO = static_cast<BS_FileOut *>(DBS->BS_FileOut);
	BS_BallScrew *BS = static_cast<BS_BallScrew *>(DBS->BS_BallScrew);

	FO->init(*FI);

	BS->allocate(*FI);
	BS->init(*FI, FI->stt.v0, FI->stt.w0, FI->stt.wn);

	Rigid::l = FI->rigid.l;
	Rigid::t = FI->rigid.t;

	DBS->n = 13 * (FI->ballnum + 2);
	DBS->y = new double[DBS->n];
	DBS->dydt = new double[DBS->n];

	return;
}

void BS_preset(DynamicBallscrew*DBS) {

	BS_BallScrew *BS = static_cast<BS_BallScrew *>(DBS->BS_BallScrew);

	BS->preset_y0(1e-9, 1e-12);

	return;
}

int BS_size_y(DynamicBallscrew*DBS) {

	BS_FileIn *FI = static_cast<BS_FileIn *>(DBS->BS_FileIn);

	return 13 * (FI->ballnum + 2);
}

void BS_set_y(DynamicBallscrew*DBS, const double*y) {

	BS_BallScrew *BS = static_cast<BS_BallScrew *>(DBS->BS_BallScrew);

	BS->set_y(y);

	return;
}

void BS_get_y(DynamicBallscrew*DBS, double*y) {

	BS_BallScrew *BS = static_cast<BS_BallScrew *>(DBS->BS_BallScrew);

	BS->get_y(y);

	return;
}

void BS_get_dydt(DynamicBallscrew*DBS, double*dydt) {

	BS_BallScrew *BS = static_cast<BS_BallScrew *>(DBS->BS_BallScrew);

	BS->get_dydt(DBS->dydt, 0.0, 0.0);

	for (int i = 0; i < DBS->n; i++)
		dydt[i] = DBS->dydt[i];

	return;
}

void deleteBS(DynamicBallscrew*DBS) {
	delete DBS;
}




