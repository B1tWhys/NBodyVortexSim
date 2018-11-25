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
 NOTE:
 GS <=> group seperator <=> 0x1D
 RS <=> record seperator <=> 0x1E

 File Format:

 Each timestep starts with a GS, followed by either "RK" or "TS" to distinguish between
 TimeSteps and RungeKutta steps.

 After the metadata line is the vortex data for that complete timestep. Each data section 
 starts with a RS. Each following line (carriage return delimited) contains the data for a
 vortex (details below). If the file contains runge-kutta steps, 4 RK steps follow the timestep line.

 Vortex data is succeeded by tracer data, which is in a similar format.

 Tracer data is followed optionally by 4 runge-kutta steps. Each RK step starts with RK, and is followed by
 a 2D array of data in python syntax. This means that the data can be parsed into arrays just by running the
 python exec() function on the array. Each 3D array is an array of 4 time steps, each timestep consisting of
 a dictionary where the keys are vIDs, and the values are arrays of x-y positions. Maybe this will make
 things more clear:

[
    # RK step 1
    {
        1: [xPosition, yPosition], # data for vortID = 1
        2: [xPosition, yPosition], # data for vortID = 2
        ...
    },
    # RK step 2
    {
        1: [xPosition, yPosition], # data for vortID = 1
        2: [xPosition, yPosition], # data for vortID = 2
        ...
    },
    ... # RK step 3 & 4
]

 
 Details:

 <GS>ST,Step #,Time,Seed Val,#vorts,#tracers
 <RS>vID,xPos,yPos,xVel,yVel,Intensity,spawnStep
 vID,xPos,yPos,xVel,yVel,Intensity,spawnStep
 vID,xPos,yPos,xVel,yVel,Intensity,spawnStep
 vID,xPos,yPos,xVel,yVel,Intensity,spawnStep
 .
 .
 .
 <RS>RK,1
 vID,xPos,yPos,xVel,yVel,Intensity,spawnStep
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
				"%li,%.15f,%.15f,%.15f,%.15f,%.15f,%i\n",
				vort->vID,
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

void saveIntermediateVortPositions(int numVorts, struct RKPositions *positions) {
    fprintf(file, "[");
    for (int RKStep = 1; RKStep <= 4; ++RKStep) {
        fprintf(file, "{");
        for (int i = 0; i < numVorts; ++i) {
            double x, y;
            switch(RKStep) {
                case 1:
                    x = positions->step1Pos.x;
                    y = positions->step1Pos.y;
                    break;
                case 2:
                    x = positions->step2Pos.x;
                    y = positions->step2Pos.y;
                    break;
                case 3:
                    x = positions->step3Pos.x;
                    y = positions->step3Pos.y;
                    break;
                case 4:
                    x = positions->step4Pos.x;
                    y = positions->step4Pos.y;
                    break;
            }
            fprintf(file, "%ld: [%lf, %lf],", positions->vort->vID, x, y);
        }
        fprintf(file, "}");
    }
    fprintf(file, "]");
}

void saveState_binary(int timestep, double currentTime, unsigned int currentSeed, int numVorts, int numTracers, struct Vortex *vorts, struct Tracer *tracers) {
	assert(fprintf(file, "\x1D%i,%f,%u,%i,%i\n", timestep, currentTime, currentSeed, numVorts, numTracers) >= 0);
}

void closeFile() {
	fprintf(stderr, "closing file\n");
	fclose(file);
}

/**
 binary search for the beginning of the target timestep
 */
long findTimeStep(FILE *f, int target, long lowPos, long highPos) {
	printf("searching in %li-%li\n", lowPos, highPos);
	char strbuff[100];
	long domainWidth = highPos - lowPos;
	fseek(f, lowPos + domainWidth/2, SEEK_SET);
	
	while (strbuff[0] != 0x1D) {
		*strbuff = (char)fgetc(f);
		// if our entire domain is *inside* of the last timestep (and doesn't include that timestep's record sep.) then
		// we will run out of file b4 we find the TS header. If this happens, then we should seek back in the file by the
		// width of our search domain, then search from there. This searches the other half of the parent's domain.
		if (*strbuff == -1) {
			if (!fseek(f, lowPos - highPos, SEEK_CUR)) {
				exit(EXIT_FAILURE); // seek failed.
			}
		}
	}
	*strbuff = 0;
	long tsStartLoc = ftell(f);
	
	for (int i = 0; 1; i++) {
		char nextChar = fgetc(f);
		if (nextChar != ',' && nextChar != '\n') {
			strbuff[i] = nextChar;
		} else {
			break;
		}
	}
	
	int tsNum = atoi(strbuff);
//	printf("center search val: %s\n", strbuff);

	if (tsNum == target) {
		return tsStartLoc;
	} else if (tsNum > target) {
		return findTimeStep(f, target, lowPos, highPos - domainWidth/2);
	} else {
		return findTimeStep(f, target, lowPos + domainWidth/2, highPos);
	}
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
	
	fseek(sourceF, 0L, SEEK_END);
	long fLen = ftell(sourceF);
	fseek(sourceF, fLen/2, SEEK_SET);
	long timestepStart = findTimeStep(sourceF, loadIndex, 0L, fLen);
	fseek(sourceF, timestepStart, SEEK_SET);
	
	readNextCSV;
	assert(loadIndex == atoi(strbuff)); // check that we are at the correct timestep in the file
	currentTimestep = loadIndex;
	
	// handle (or don't handle) time
//	clearStrBuff;
	readNextCSV;
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
		vort->vID = atoi(strbuff);
        vort->vIndex = vIndex;
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
	fclose(sourceF);
	printf("file read finished\n");
}
