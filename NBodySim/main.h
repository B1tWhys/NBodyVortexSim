//
//  main.h
//  NBodySim
//
//

#ifndef main_h
#define main_h

#include <pthread.h>

struct Vortex {
	int vIndex;
	int initStep;
	
	double intensity;
	double *position;
	double *velocity; // this is change in coord. per time step
	volatile int otherVort;
};

struct Tracer {
	int tIndex;
	
	double *position;
	double *velocity; // this is change in coord. per time step
};

struct TracerArgs {
	int RKStep;
	struct Tracer *tracers;
	int numTracers;
	double *tracerRadii;
	double *intermediateTracerRads;
	struct Vortex *vortices;
};

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

#endif /* main_h */
