//
//  main.h
//  NBodySim
//
//  Created by Skyler Arnold on 6/25/18.
//  Copyright © 2018 Skyler Arnold. All rights reserved.
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

#endif /* main_h */
