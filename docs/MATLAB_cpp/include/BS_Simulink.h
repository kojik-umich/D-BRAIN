#ifndef BS_SIMULINK_H
#define BS_SIMULINK_H

#ifdef _WIN32
#ifdef EXPORT
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT __declspec(dllimport)
#endif
#else
#define DLL_EXPORT
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

DLL_EXPORT struct BS_Simulink* new_BS(void);
DLL_EXPORT void BS_initialize(struct BS_Simulink* Simulink, int n);
DLL_EXPORT int BS_sum(struct BS_Simulink*Simulink);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif

