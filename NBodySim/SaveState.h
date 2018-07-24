//
//  SaveState.h
//  NBodySim
//
//  Created by Skyler Arnold on 7/23/18.
//

#ifndef SaveState_h
#define SaveState_h

#include "main.h"

#include <stdio.h>

void saveState(int timestep, double currentTime, unsigned int currentSeed, int numVorts, int numTracers, struct Vortex *vorts, struct Tracer *tracers);

#endif /* SaveState_h */
