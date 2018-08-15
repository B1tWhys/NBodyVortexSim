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
#include <stdio.h>

void openFile(void);
void saveState(int timestep, double currentTime, long currentSeed, int numVorts, int numTracers, struct Vortex *vorts, struct Tracer *tracers);
void closeFile(void);

#endif /* SaveState_h */
