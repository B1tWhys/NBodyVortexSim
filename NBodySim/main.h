//
//  main.h
//  NBodySim
//
//

#ifndef main_h
#define main_h

#include <pthread.h>

extern int currentTimestep;

struct Vortex {
    long vID; // unique vortex ID
	int vIndex;
	int initStep;
	
	double intensity;
	double *position;
	double *velocity; // this is change in coord. per time step
};

struct Vector {
    double x;
    double y;
    double mag;
};

struct RKPositions {
    struct Vortex *vort;
    struct Vector step1Pos;
    struct Vector step2Pos;
    struct Vector step3Pos;
    struct Vector step4Pos;
};

struct Tracer {
	int tIndex;
	
	double *position;
	double *velocity; // this is change in coord. per time step
};

// used to pass vortex info between threads
struct VortexArgs {
	int RKStep;
	struct Vortex *vortices;
	int originVortIndex;
	double *intermediateRadii;
	double *workingRadii;
	long vortRadLen;
	double *intermediateTracerRads;
	int numTracers;
};

// used to pass tracer info between threads
struct TracerArgs {
	int RKStep;
	struct Tracer *tracers;
	int numTracers;
	double *tracerRadii;
	double *intermediateTracerRads;
	struct Vortex *vortices;
};

#endif /* main_h */
