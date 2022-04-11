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
	void *BS;
};
typedef struct BS_Simulink Simulink;

DLL_EXPORT BS_Simulink*newSimulink();
DLL_EXPORT void BS_initialize(BS_Simulink* Simulink, int n);
DLL_EXPORT int BS_sum(BS_Simulink*Simulink);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif

