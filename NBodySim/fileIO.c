//
//  SaveState.c
//  NBodySim
//

#include "fileIO.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#ifndef DATA_OUT_FILEPATH
#define DATA_OUT_FILEPATH "./data/rawData"
#endif

/*
 File Format:
 
 <GS>Step #,Time,Seed Val,#vorts,#tracers
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
	// file = fopen("./data/rawData", "a");
	file = fopen(DATA_OUT_FILEPATH, "w");
	assert(file);
}

void saveState(int timestep, long currentSeed, int numVorts, int numTracers, struct Vortex *vorts, struct Tracer *tracers) {
	assert(fprintf(file, "\x1D%i,%li,%i,%i\n", timestep, currentSeed, numVorts, numTracers) >= 0);
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
	assert(fprintf(file, "\x1D%i,%f,%u,%i,%i\n", timestep, currentTime, currentSeed, numVorts, numTracers) >= 0);
}

void closeFile() {
	fprintf(stderr, "closing file\n");
	fclose(file);
}

// write 0's into strbuff. This assumes strbuff is 100 chars long.
#define clearStrBuff for (int i = 0; i < 100; i++) strbuff[i] = 0

// read characters into strbuff until we get to the next ',' or '\n'
#define readNextCSV \
for (int i = 0; 1; i++) {\
char nextChar = fgetc(sourceF);\
if (nextChar != ',' && nextChar != '\n') {\
strbuff[i] = nextChar;\
} else {\
break;\
}\
}

void initFromFile(char *fName, int loadIndex, struct Vortex *vortices[], int *numDriverVorts, int *vortsAllocated, struct Tracer *tracers[]) {
	FILE *sourceF = fopen(fName, "r");
	char strbuff[100]; // if a number in the file exceeds 100 characters in length, this will overflow
	clearStrBuff;
	
	// loop through timesteps until we are at the timestep to load from
	for (int tsIndex = 0; tsIndex <= loadIndex; tsIndex++) {
		while (strbuff[0] != 0x1D) { // loop until we reach the beginning of the next timestep
			strbuff[0] = (char)fgetc(sourceF);
			assert(strbuff[0] != -1); // we reached EOF before reaching the designated timestep
		}
		strbuff[0] = 0;
	}
	
	readNextCSV;
	assert(loadIndex == atoi(strbuff)); // check that we are at the correct timestep in the file
	currentTimestep = loadIndex;
	
	// handle (or don't handle) time
//	clearStrBuff;
//	readNextCSV;
//	*time = atof(strbuff);
	
	clearStrBuff;
	readNextCSV;
	lastX = atol(strbuff);
	
	clearStrBuff;
	readNextCSV;
	*numDriverVorts = atoi(strbuff);
	
	clearStrBuff;
	readNextCSV;
	assert(NUM_TRACERS == atoi(strbuff)); // the number of tracers in the file != num tracers in the current config file
	
	clearStrBuff;
	strbuff[0] = fgetc(sourceF);
	assert(strbuff[0] == 0x1E);
	
	*vortsAllocated = (*numDriverVorts) * 1.5;
	*vortices = malloc((int)sizeof(struct Vortex) * (*vortsAllocated));
	*tracers = malloc((int)sizeof(struct Vortex) * NUM_TRACERS * 1.5);
	
	for (int vIndex = 0; vIndex < *numDriverVorts; vIndex++) {
		struct Vortex *vort = &(*vortices)[vIndex];
		vort->position = malloc(sizeof(double) * 2);
		vort->velocity = malloc(sizeof(double) * 2);
		clearStrBuff;
		readNextCSV;
		vort->vIndex = atoi(strbuff);
		clearStrBuff;
		readNextCSV;
		vort->position[0] = atof(strbuff);
		clearStrBuff;
		readNextCSV;
		vort->position[1] = atof(strbuff);
		clearStrBuff;
		readNextCSV;
		vort->velocity[0] = atof(strbuff);
		clearStrBuff;
		readNextCSV;
		vort->velocity[1] = atof(strbuff);
//		readNextCSV; // skip over totVel
		clearStrBuff;
		readNextCSV;
		vort->intensity = atof(strbuff);
		clearStrBuff;
		readNextCSV;
		vort->initStep = atoi(strbuff);
	}
	
	assert(fgetc(sourceF) == 0x1E); // make sure we are past the vortices
	
	for (int tIndex = 0; tIndex < NUM_TRACERS; tIndex++) {
		struct Tracer *tracer = &(*tracers)[tIndex];
		tracer->position = malloc(sizeof(double) * 2);
		tracer->velocity = malloc(sizeof(double) * 2);
		clearStrBuff;
		readNextCSV;
		tracer->tIndex = atoi(strbuff);
		clearStrBuff;
		readNextCSV;
		tracer->position[0] = atof(strbuff);
		clearStrBuff;
		readNextCSV;
		tracer->position[1] = atof(strbuff);
		clearStrBuff;
		readNextCSV;
		tracer->velocity[0] = atof(strbuff);
		clearStrBuff;
		readNextCSV;
		tracer->velocity[1] = atof(strbuff);
//		readNextCSV; // skip totVel
	}
	
	assert(fgetc(sourceF) == 0x1D); // make sure that's the last of the tracers for that timestep.
}
