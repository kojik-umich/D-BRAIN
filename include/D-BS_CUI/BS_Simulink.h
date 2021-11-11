#ifndef BS_WRAPPER_H
#define BS_WRAPPER_H

#ifdef EXPORT
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT __declspec(dllimport)
#endif

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

struct BS_Simulink {
	void *Ball0;
	void *Ball1;
	void *BS_Nut;
	void *BS_Shaft;
	void *BS_BallNutPair0;
	void *BS_BallNutPair1;
	void *BS_BallShaftPair0;
	void *BS_BallShaftPair1;
	void *BS_SimulinkIn;
};
typedef struct BS_Simulink Simulink;

DLL_EXPORT BS_Simulink * newSimulink();
DLL_EXPORT void BS_read(BS_Simulink*Simulink, const char dbsin_csv[]);
DLL_EXPORT void BS_init(BS_Simulink*Simulink);
DLL_EXPORT void deleteBS(BS_Simulink*Simulink);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif

