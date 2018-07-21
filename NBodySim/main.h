//
//  main.h
//  NBodySim
//
//

#ifndef main_h
#define main_h

// A vortex is either active, meaning that it sould continue existing
// or its state is Merged meaning that it was merged into another vortex, and
// is awaiting being either deleted or re-initialized at a new location with
// a new intensity
//enum vortStatus {Active, Merged};

struct Vortex {
	int vIndex;
//	enum vortStatus status;
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
