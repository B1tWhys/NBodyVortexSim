//
//  SaveState.c
//  NBodySim
//

#include "SaveState.h"
#include <stdio.h>
#include <assert.h>

/*
 File Format:
 
 <GS>Step #,Time,Seed Val,#vorts,#tracers
 <RS>vIndex,xPos,yPos,xVel,yVel,totVel,Intensity,spawnStep
 vIndex,xPos,yPos,xVel,yVel,totVel,Intensity,spawnStep
 vIndex,xPos,yPos,xVel,yVel,totVel,Intensity,spawnStep
 vIndex,xPos,yPos,xVel,yVel,totVel,Intensity,spawnStep
 .
 .
 .
 <RS>tIndex,xPos,yPos,xVel,yVel,totVel
 tIndex,xPos,yPos,xVel,yVel,totVel
 tIndex,xPos,yPos,xVel,yVel,totVel
 tIndex,xPos,yPos,xVel,yVel,totVel
 .
 .
 .
 <GS>Step #,Time,Seed Val,#vorts,#tracers
 .
 .
 .
 */

FILE *file;

void openFile() {
	// file = fopen("./data/rawData", "a");
	file = fopen(DATA_FILEPATH, "w");
	assert(file);
}

void saveState(int timestep, double currentTime, unsigned int currentSeed, int numVorts, int numTracers, struct Vortex *vorts, struct Tracer *tracers) {
	assert(fprintf(file, "\x1D%i,%f,%u,%i,%i\n", timestep, currentTime, currentSeed, numVorts, numTracers) >= 0);
	fputc(0x1E, file);
	for (int i = 0; i < numVorts; i++) {
		struct Vortex *vort = &vorts[i];
		assert(fprintf(file,
				"%i,%.15f,%.15f,%.15f,%.15f,%.15f,%i\n",
				vort->vIndex,
				vort->position[0],
				vort->position[1],
				vort->velocity[0],
				vort->velocity[1],
				vort->intensity,
				vort->initStep) >= 0);
	}
	fputc(0x1E, file);
	
	for (int i = 0; i < numTracers; i++) {
		struct Tracer *tracer = &tracers[i];
		assert(fprintf(file,
				"%i,%.15f,%.15f,%.15f,%.15f\n",
				tracer->tIndex,
				tracer->position[0],
				tracer->position[1],
				tracer->velocity[0],
				tracer->velocity[1]) >= 0);
	}
	
	// fputc(29, file);
	fflush(file);
}

void saveState_binary(int timestep, double currentTime, unsigned int currentSeed, int numVorts, int numTracers, struct Vortex *vorts, struct Tracer *tracers) {
	void *writePtr = &timestep;
	assert(fprintf(file, "\x1D%i,%f,%u,%i,%i\n", timestep, currentTime, currentSeed, numVorts, numTracers) >= 0);
	
}

void closeFile() {
	fprintf(stderr, "closing file\n");
	fclose(file);
}
