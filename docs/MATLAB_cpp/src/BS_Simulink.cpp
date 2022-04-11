#define EXPORT

#include "BS.h"
#include "BS_Simulink.h"

DLL_EXPORT BS_Simulink *newSimulink() {

	BS_Simulink*Simulink = new BS_Simulink();

	Simulink->BS = new BS();

	return Simulink;
}

DLL_EXPORT void BS_initialize(BS_Simulink* Simulink, int n)
{
	BS* BS_ = static_cast<BS*>(Simulink->BS);

	BS_->initialize(n);

	return;
}

DLL_EXPORT int BS_sum(BS_Simulink* Simulink)
{
	BS* BS_ = static_cast<BS*>(Simulink->BS);

	return BS_->sum();
}

DLL_EXPORT void deleteBS(BS_Simulink*Simulink) {
	delete Simulink;
}




