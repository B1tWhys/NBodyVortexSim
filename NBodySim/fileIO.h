//
//  SaveState.h
//  NBodySim
//
//  Created by Skyler Arnold on 7/23/18.
//

#ifndef SaveState_h
#define SaveState_h

#include "main.h"
#include "constants.h"
#include "RNG.h"
#include <stdio.h>

void openFile(void);
void saveState(int timestep, long currentSeed, int numVorts, int numTracers, struct Vortex *vorts, struct Tracer *tracers);
void saveRKState(int RKStep, int timestep, long currentSeed, int numVorts, int numTracers, struct Vortex *vorts, struct Tracer *tracers);
void closeFile(void);

void initFromFile(char *fName, int loadIndex, struct Vortex *vortices[], int *numDriverVorts, int *vortsAllocated, struct Tracer *tracers[]);

#endif /* SaveState_h */
