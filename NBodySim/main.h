//
//  main.h
//  NBodySim
//
//

#ifndef main_h
#define main_h

struct Vortex {
	int vIndex;
	int initStep;
	
	double intensity;
	double *position;
	double *velocity; // this is change in coord. per time step
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

#endif /* main_h */
