//
//  SaveState.c
//  NBodySim
//

#include "SaveState.h"

#include <stdio.h>

/*
 File Format:
 
 Step #,Time,Seed Val,#vorts,#tracers
 <RS>vIndex,xPos,yPos,xVel,yVel,Intensity,spawnStep
 vIndex,xPos,yPos,xVel,yVel,Intensity,spawnStep
 vIndex,xPos,yPos,xVel,yVel,Intensity,spawnStep
 vIndex,xPos,yPos,xVel,yVel,Intensity,spawnStep
 .
 .
 .
 <RS>tIndex,xPos,yPos,xVel,yVel
 tIndex,xPos,yPos,xVel,yVel
 tIndex,xPos,yPos,xVel,yVel
 tIndex,xPos,yPos,xVel,yVel
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
	file = fopen("./data/rawData", "a");
}

void saveState(int timestep, double currentTime, unsigned int currentSeed, int numVorts, int numTracers, struct Vortex *vorts, struct Tracer *tracers) {
	fprintf(file, "%i,%f,%u,%i,%i\n", timestep, currentTime, currentSeed, numVorts, numTracers);
	fputc(30, file);
	for (int i = 0; i < numVorts; i++) {
		struct Vortex *vort = &vorts[i];
		fprintf(file,
				"%i,%.15f,%.15f,%.15f,%.15f,%.15f,%i\n",
				vort->vIndex,
				vort->position[0],
				vort->position[1],
				vort->velocity[0],
				vort->velocity[1],
				vort->intensity,
				vort->initStep);
	}
	fputc(30, file);
	
	for (int i = 0; i < numTracers; i++) {
		struct Tracer *tracer = &tracers[i];
		fprintf(file,
				"%i,%.15f,%.15f,%.15f,%.15f\n",
				tracer->tIndex,
				tracer->position[0],
				tracer->position[1],
				tracer->velocity[0],
				tracer->velocity[1]);
	}
	
	fputc(29, file);
}

void closeFile() {
	fclose(file);
}
