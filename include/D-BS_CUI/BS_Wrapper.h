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

struct DynamicBallscrew {
	void *BS_BallScrew;
	void *BS_FileIn;
	void *BS_FileOut;
	double *y;
	double *dydt;
	int n;
};
typedef struct DynamicBallscrew DBS;

DLL_EXPORT DynamicBallscrew * newBS();
DLL_EXPORT void BS_read(DynamicBallscrew*DBS, const char dbsin_csv[]);
DLL_EXPORT void BS_init(DynamicBallscrew*DBS);
DLL_EXPORT void BS_preset(DynamicBallscrew*DBS);
DLL_EXPORT int BS_size_y(DynamicBallscrew*DBS);
DLL_EXPORT void BS_set_y(DynamicBallscrew*DBS, const double*y);
DLL_EXPORT void BS_get_y(DynamicBallscrew*DBS, double*y);
DLL_EXPORT void BS_get_dydt(DynamicBallscrew*DBS, double*dydt);
DLL_EXPORT void deleteBS(DynamicBallscrew*DBS);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif

