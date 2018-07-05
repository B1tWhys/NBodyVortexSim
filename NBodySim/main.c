//
//  main.c
//  NBodySim
//
//  Created by Skyler Arnold on 6/12/18.
//  Copyright © 2018 Skyler Arnold. All rights reserved.
//

#include "main.h"
#include "guiOutput.h"
#include "constants.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <string.h>

unsigned int randomSeed;
int currentTimestep;

#pragma mark - Index Calculators

long calculateVortexRadiiIndex(long vortIndex1, long vortIndex2) {
	// radii is basically a hash table, with this as the hash function.
	// I suspect that this is better than the standard dictionary for this.
	// I'm able to guarentee that the structure is full, and that there are no
	// collisions, and I'm able to iterate efficiently by incrementing an array pointer.
	if (vortIndex1 < vortIndex2) {
		return ((vortIndex2-1)*vortIndex2/2 + vortIndex1) * 3;
	} else {
		return ((vortIndex1-1)*vortIndex1/2 + vortIndex2) * 3;
	}
}

long calculateTracerRadiiIndex(long tracerIndex, long vortIndex) {
	/*
	 rows = tracers
	 cols = vorts
	 */

	return (vortIndex * NUM_TRACERS + tracerIndex)*3;
}

#pragma mark - Diff Eq code

// this function doesn't work, so I'm using the pythagorean solver
void updateRadii_derivative(double *radii, struct Vortex *vortices, int numvortices) {
	long index;

	for (int i = 1; i < numvortices; ++i) {
		index = i*(i+1)/2;
		for (int j = 0; j < i; ++j) {
			radii[index] += (radii[index+1]*(vortices[i].velocity[1] - vortices[j].velocity[1]) +
							 radii[index+2]*(vortices[i].velocity[0] - vortices[j].velocity[0])) / radii[index];
			radii[index+1] = vortices[i].position[0] - vortices[j].position[0];
			radii[index+2] = vortices[i].position[1] - vortices[j].position[1];
			++index;
		}
	}
}

// This works reliably but it *may* be *very slightly* slower, and may also be completely redundant b/c of the RK4 function
void updateRadii_pythagorean(double *vortexRadii, struct Vortex *vortices, int numVortices, double *tracerRadii, struct Tracer *tracers, int numTracers) {
	long index;

	for (int i = 0; i < numVortices; ++i) { // if this doesn't segfault it'll be a goddamn miracle
		for (int j = 0; j < i; ++j) {
			index = calculateVortexRadiiIndex(vortices[i].vIndex, vortices[j].vIndex);

			vortexRadii[index+1] = vortices[i].position[0] - vortices[j].position[0];
			vortexRadii[index+2] = vortices[i].position[1] - vortices[j].position[1];
			vortexRadii[index] = sqrt(pow(vortexRadii[index+1], 2) + pow(vortexRadii[index+2], 2));
		}
	}
	index = 0;
	for (int tracerIndex = 0; tracerIndex < numTracers; tracerIndex++) {
		struct Tracer *tracer = &tracers[tracerIndex];
		for (int vortIndex = 0; vortIndex < numVortices; vortIndex++) {
			struct Vortex *vort = &vortices[vortIndex];
			index = calculateTracerRadiiIndex(tracerIndex, vortIndex);

			tracerRadii[index+1] = vort->position[0] - tracer->position[0];
			tracerRadii[index+2] = vort->position[1] - tracer->position[1];
			tracerRadii[index] = sqrt(pow(tracerRadii[index+1], 2) + pow(tracerRadii[index+2], 2));
		}
	}
}

double velocityFunc(double vortex2Intensity, double radius) { // calculates velocity for vortex 1 as a result of vortex 2
	return vortex2Intensity/(2.*M_PI*radius);
}

char domains = 0;

void calculateDxDj_DyDj_vortex(double *dxdj, double *dydj, struct Vortex *vort, struct Vortex *vortices, int numVorts, double *intRads, long numIntRads) {
	for (int j = 0; j < numVorts; ++j) {
		if (vort->vIndex == j) {
			continue;
		}

		long radiiIndex = calculateVortexRadiiIndex(vort->vIndex, vortices[j].vIndex);
		long intRadiiIndex = radiiIndex;

		for (int domain = 0; domain <= domains; domain++) {
			double rad;
			double xRad = (vort->vIndex < vortices[j].vIndex) ? intRads[intRadiiIndex + 1] : -intRads[intRadiiIndex + 1];
			double yRad = (vort->vIndex < vortices[j].vIndex) ? intRads[intRadiiIndex + 2] : -intRads[intRadiiIndex + 2];

			if (domain) {
				switch (domain) {
					case 1: {
						xRad -= DOMAIN_SIZE_X;
						yRad += DOMAIN_SIZE_Y;
						break;
					} case 2: {
						yRad += DOMAIN_SIZE_Y;
						break;
					} case 3: {
						xRad += DOMAIN_SIZE_X;
						yRad += DOMAIN_SIZE_Y;
						break;
					} case 4: {
						xRad -= DOMAIN_SIZE_X;
						break;
					} case 5: {
						xRad += DOMAIN_SIZE_X;
						break;
					} case 6: {
						xRad -= DOMAIN_SIZE_X;
						yRad -= DOMAIN_SIZE_Y;
						break;
					} case 7: {
						yRad -= DOMAIN_SIZE_Y;
						break;
					} case 8: {
						xRad += DOMAIN_SIZE_X;
						yRad -= DOMAIN_SIZE_Y;
						break;
					} default: {
						break;
					}
				}
				rad = sqrt(pow(xRad, 2) + pow(yRad, 2));
			} else {
				rad = intRads[intRadiiIndex];
			}

			if (rad > DOMAIN_SIZE_X) {
				continue; // domain truncation
			}

			double vmag = velocityFunc(vortices[j].intensity, rad); // total velocity magnitude
			*dxdj +=  (yRad/rad) * TIMESTEP * vmag;
			*dydj += (-xRad/rad) * TIMESTEP * vmag;
		}
	}
}

void calculateDxDj_DyDj_tracer(double *dxdj, double *dydj, struct Tracer *tracer, struct Tracer *tracers, int numTravers, double *intRads, long numIntRads, struct Vortex *vortices, int numVorts) {
	for (int vortIndex = 0; vortIndex < numVorts; vortIndex++) {
		long intRadIndex = calculateTracerRadiiIndex(tracer->tIndex, vortIndex);

		for (int domain = 0; domain <= domains; domain++) {
			double rad;
			double xRad = intRads[intRadIndex + 1], yRad = intRads[intRadIndex + 2];

			if (domain) {
				switch (domain) {
					case 1: {
						xRad -= DOMAIN_SIZE_X;
						yRad += DOMAIN_SIZE_Y;
						break;
					} case 2: {
						yRad += DOMAIN_SIZE_Y;
						break;
					} case 3: {
						xRad += DOMAIN_SIZE_X;
						yRad += DOMAIN_SIZE_Y;
						break;
					} case 4: {
						xRad -= DOMAIN_SIZE_X;
						break;
					} case 5: {
						xRad += DOMAIN_SIZE_X;
						break;
					} case 6: {
						xRad -= DOMAIN_SIZE_X;
						yRad -= DOMAIN_SIZE_Y;
						break;
					} case 7: {
						yRad -= DOMAIN_SIZE_Y;
						break;
					} case 8: {
						xRad += DOMAIN_SIZE_X;
						yRad -= DOMAIN_SIZE_Y;
						break;
					} default: {
						break;
					}
				}
				rad = sqrt(pow(xRad, 2) + pow(yRad, 2));
			} else {
				rad = intRads[intRadIndex];
			}
			
			if (rad > DOMAIN_SIZE_X) {
				continue;
			}

			double vmag = velocityFunc(vortices[vortIndex].intensity, rad);
			double dxInc = (yRad/rad) * TIMESTEP * vmag;
			double dyInc = (-xRad/rad) * TIMESTEP * vmag;
			*dxdj +=  (yRad/rad) * TIMESTEP * vmag;
			*dydj += (-xRad/rad) * TIMESTEP * vmag;
		}
	}
}

void stepForward_RK4(struct Vortex *vortices, double *vortRadii, int numVortices, double *tracerRadii, struct Tracer *tracers, int numTracers) {
	/*
	 There is definately some duplicate computation being done here.     ୧(ಠ Д ಠ)୨

	  this algorithm ends up making 2 coppies of each radius. One in each direction.

	 // upon further consideration this was very stupid. Definately just do the whole matrix

	 *also im not sure if i am allowed to do this summing of k terms...?
	 */

	int sizeOfRadEntry = sizeof(double) * 3;
	long intermediateVortRadSize = (pow(numVortices, 2)-numVortices)/2 * sizeOfRadEntry;
	double *intermediateVortRads = malloc(intermediateVortRadSize);
	memcpy(intermediateVortRads, vortRadii, intermediateVortRadSize);

	long intermediateTracerRadSize = numVortices * numTracers * sizeOfRadEntry;
	double *intermediateTracerRads = malloc(intermediateTracerRadSize);
	memcpy(intermediateTracerRads, tracerRadii, intermediateTracerRadSize);

	for (int i = 0; i < numVortices; i++) {
		struct Vortex *vort = &vortices[i];
		vort->velocity[0] = 0;
		vort->velocity[1] = 0;
	}

	for (int i = 0; i < numTracers; i++) {
		struct Tracer *tracer = &tracers[i];
		tracer->velocity[0] = 0;
		tracer->velocity[1] = 0;
	}

	for (int RKStep = 1; RKStep <= 4; RKStep++) {
		double dxdj = 0, dydj = 0;

		for (int originVortIndex = 0; originVortIndex < numVortices; ++originVortIndex) {
			double k1_x = 0, k2_x = 0, k3_x = 0, k4_x = 0;
			double k1_y = 0, k2_y = 0, k3_y = 0, k4_y = 0;

			struct Vortex *vort = &vortices[originVortIndex];
			int vortIndex = vort->vIndex;

			dxdj = 0;
			dydj = 0;

			calculateDxDj_DyDj_vortex(&dxdj, &dydj, vort, vortices, numVortices, intermediateVortRads, numVortices);

			for (int j = 0; j < numVortices; ++j) {
				if (originVortIndex == j) {
					continue;
				}
				
				long radiiIndex = calculateVortexRadiiIndex(vortIndex, vortices[j].vIndex);
				long intRadiiIndex = radiiIndex;

				switch (RKStep) {
					case 1: {
						k1_x += dxdj;
						k1_y += dydj;

						dxdj = dxdj/2.;
						dydj = dydj/2.;
						break;
					}
					case 2: {
						k2_x += dxdj;
						k2_y += dydj;

						dxdj = dxdj/2.;
						dydj = dydj/2.;
						break;
					}
					case 3: {
						k3_x += dxdj;
						k3_y += dydj;

						break;
					}
					case 4: {
						k4_x += dxdj;
						k4_y += dydj;
						break;
					}
					default:
						break;
				}

				// perhaps an optimization can be made where we do some kind of lazy copy of radii into intermediateRads,
				// initially making intermediate Rads only a pointer to radii, then doing the copy here in RK step 1.
				// This would save 1 large matrix operation.

				intermediateVortRads[intRadiiIndex+1] = vortRadii[radiiIndex + 1] + dxdj;
				intermediateVortRads[intRadiiIndex+2] = vortRadii[radiiIndex + 2] + dydj;

				// pythagorean.
				intermediateVortRads[intRadiiIndex] = sqrt(pow(intermediateVortRads[intRadiiIndex+1], 2) + pow(intermediateVortRads[intRadiiIndex+2], 2));
			}
			vort->velocity[0] += (k1_x + k2_x * 2 + k3_x * 2 + k4_x)/6;
			vort->velocity[1] += (k1_y + k2_y * 2 + k3_y * 2 + k4_y)/6;
		}

		for (int tracerIndex = 0; tracerIndex < numTracers; ++tracerIndex) {
			double k1_x = 0, k2_x = 0, k3_x = 0, k4_x = 0;
			double k1_y = 0, k2_y = 0, k3_y = 0, k4_y = 0;

			struct Tracer *tracer = &tracers[tracerIndex];
			int tracerIndex = tracer->tIndex; // TODO: straight up delete this line. dont need to change anything anywhere
			
			dxdj = 0;
			dydj = 0;


			calculateDxDj_DyDj_tracer(&dxdj,
									  &dydj,
									  tracer,
									  tracers,
									  numTracers,
									  intermediateTracerRads,
									  (numVortices * numTracers),
									  vortices,
									  numVortices);

			for (int vortexIndex = 0; vortexIndex < numVortices; ++vortexIndex) {
				long radIndex = calculateTracerRadiiIndex(tracerIndex, vortexIndex);

				switch (RKStep) {
					case 1: {
						k1_x += dxdj;
						k1_y += dydj;

						dxdj = dxdj/2;
						dydj = dydj/2;
						break;
					}
					case 2: {
						k2_x += dxdj;
						k2_y += dydj;

						dxdj = dxdj/2;
						dydj = dydj/2;
						break;
					}
					case 3: {
						k3_x += dxdj;
						k3_y += dydj;

						break;
					}
					case 4: {
						k4_x += dxdj;
						k4_y += dydj;
						break;
					}
					default:
						break;
				}

				// perhaps an optimization can be made where we do some kind of lazy copy of radii into intermediateRads,
				// initially making intermediate Rads only a pointer to radii, then doing the copy here in RK step 1.
				// This would save 1 large matrix operation.

				intermediateTracerRads[radIndex + 1] = tracerRadii[radIndex + 1] + dxdj;
				intermediateTracerRads[radIndex + 2] = tracerRadii[radIndex + 2] + dydj;

				// pythagorean.
				intermediateTracerRads[radIndex] = sqrt(pow(intermediateTracerRads[radIndex+1], 2) + pow(intermediateTracerRads[radIndex+2], 2));
			}
			tracer->velocity[0] += (k1_x + k2_x * 2. + k3_x * 2. + k4_x)/6.;
			tracer->velocity[1] += (k1_y + k2_y * 2. + k3_y * 2. + k4_y)/6.;
		}
	}

	for (int i = 0; i < numVortices; ++i) {
		struct Vortex *vort = &vortices[i];
		vort->position[0] += vort->velocity[0];
		vort->position[1] += vort->velocity[1];
	}

	for (int i = 0; i < numTracers; ++i) {
		struct Tracer *tracer = &tracers[i];
		tracer->position[0] += tracer->velocity[0];
		tracer->position[1] += tracer->velocity[1];
	}

	free(intermediateVortRads);
	free(intermediateTracerRads);
}

double generateRandInRange(double lowerBound, double upperBound) { // range is inclusive on both ends
	srand(randomSeed++);
	double range = upperBound-lowerBound + 1;
	double result = (double)random()/RAND_MAX*(range+1) + lowerBound;
	return result;
}

#pragma mark - Designated Initializers

double minRad(double *radArr, long numPts) {
	if (numPts == 0) return 0;
	long radLen = (pow(numPts, 2)-numPts)/2 * 3;

	double min = radArr[0];
	for (int i = 0; i < radLen; i += 3) {
		if (radArr[i] < min) min = radArr[i];
	}
	return min;
}

void deleteVortex(struct Vortex *vort, double *vortexRadii, int *numActiveVorts, struct Vortex *vorts, double *tracerRads) {
	int deletionIndex = vort->vIndex;

	double *destPtr;
	double *sourcePtr;

	for (int rowIndex = deletionIndex; rowIndex < *numActiveVorts; rowIndex++) {
		destPtr = vortexRadii + calculateVortexRadiiIndex(0, rowIndex);
		sourcePtr = vortexRadii + calculateVortexRadiiIndex(0, rowIndex+1);

		for (int colIndex = 0; colIndex < rowIndex; colIndex++) {
			if (colIndex == deletionIndex) {
				sourcePtr += 3;
			}

			for (int i = 0; i < 3; i++) {
				*destPtr = *sourcePtr;
				destPtr++;
				sourcePtr++;
			}
		}
	}

	free(vort->position);
	free(vort->velocity);
	memcpy(&vorts[deletionIndex], &vorts[deletionIndex+1], (int)sizeof(struct Vortex) * (*numActiveVorts-deletionIndex-1));
	(*numActiveVorts)--;
	for (int vortIndex = 0; vortIndex < *numActiveVorts; vortIndex++) {
		vorts[vortIndex].vIndex--;
	}

	// TODO: remove vortex from tracer radii array


}

void mergeVorts(double *vortexRadii, int *numActiveVorts, struct Vortex *vorts, double *tracerRads) {
	for (int vort1Index = 1; vort1Index < *numActiveVorts; vort1Index++) {
		for (int vort2Index = 0; vort2Index < vort1Index; vort2Index++) {
			long radIndex = calculateVortexRadiiIndex(vort1Index, vort2Index);
			if (vortexRadii[radIndex] < VORTEX_MERGE_RADIUS_CUTOFF) {
				deleteVortex(&vorts[vort2Index], vortexRadii, numActiveVorts, vorts, tracerRads);
			}
		}
	}
}

#pragma mark - initializers

// initialize the driver vortices (in the 1st domain) with unique VIDs and random strengths and positions
void initialize_drivers_random(struct Vortex *vortices, int n, int startingID) {
	for (int i = 0; i < n; ++i) {
		vortices[i].vIndex = startingID++;
		vortices[i].position = malloc(sizeof(double) * 2);
		vortices[i].position[0] = generateRandInRange(0, DOMAIN_SIZE_X);
		vortices[i].position[1] = generateRandInRange(0, DOMAIN_SIZE_Y);
		vortices[i].velocity = calloc(2, sizeof(double));
		vortices[i].intensity = generateRandInRange(-VORTEX_INTENSITY_INIT_UPPER_BOUND, VORTEX_INTENSITY_INIT_UPPER_BOUND);
		vortices[i].initStep = 0;
	}
}

void initialize_pair_orbit_test(struct Vortex *vortices, int n) {
	assert(NUM_VORT_INIT == 2);
	assert(n == NUM_VORT_INIT);
	int intensity = VORTEX_INTENSITY_INIT_UPPER_BOUND;

	vortices[0].vIndex = 0;
	vortices[0].velocity = malloc(sizeof(double) * 2);
	vortices[0].position = malloc(sizeof(double) * 2);
	vortices[0].position[0] = DOMAIN_SIZE_X/2 - 4;
	vortices[0].position[1] = DOMAIN_SIZE_Y/2;
	vortices[0].intensity = intensity;
	vortices[0].initStep = 0;

	vortices[1].vIndex = 1;
	vortices[1].velocity = malloc(sizeof(double) * 2);
	vortices[1].position = malloc(sizeof(double) * 2);
	vortices[1].position[0] = DOMAIN_SIZE_X/2 + 4;
	vortices[1].position[1] = DOMAIN_SIZE_Y/2;
	vortices[1].intensity = intensity;
	vortices[1].initStep = 0;
}

void initialize_pair_parallel_test(struct Vortex *vortices, int n) {
	assert(NUM_VORT_INIT == 2);
	assert(n == NUM_VORT_INIT);

	double vortMagnitude = VORTEX_INTENSITY_INIT_UPPER_BOUND;

	vortices[0].vIndex = 0;
	vortices[0].velocity = malloc(sizeof(double) * 2);
	vortices[0].position = malloc(sizeof(double) * 2);
	vortices[0].position[0] = DOMAIN_SIZE_X/2 - 2;
	vortices[0].position[1] = DOMAIN_SIZE_Y/2;
	vortices[0].intensity = -vortMagnitude;
	vortices[0].initStep = 0;

	vortices[1].vIndex = 1;
	vortices[1].velocity = malloc(sizeof(double) * 2);
	vortices[1].position = malloc(sizeof(double) * 2);
	vortices[1].position[0] = DOMAIN_SIZE_X/2 + 2;
	vortices[1].position[1] = DOMAIN_SIZE_Y/2;
	vortices[1].intensity = vortMagnitude;
	vortices[1].initStep = 0;
}

void initialize_square_system(struct Vortex *vortices, int n) {
	assert(NUM_VORT_INIT == 4);
	assert(n == NUM_VORT_INIT);

	double vortMagnitude = VORTEX_INTENSITY_INIT_UPPER_BOUND;
	double offset = 10;

	vortices[0].vIndex = 0;
	vortices[0].velocity = malloc(sizeof(double) * 2);
	vortices[0].position = malloc(sizeof(double) * 2);
	vortices[0].position[0] = DOMAIN_SIZE_X/2 - offset;
	vortices[0].position[1] = DOMAIN_SIZE_Y/2 + offset;
	vortices[0].intensity = vortMagnitude;
	vortices[0].initStep = 0;

	vortices[1].vIndex = 1;
	vortices[1].velocity = malloc(sizeof(double) * 2);
	vortices[1].position = malloc(sizeof(double) * 2);
	vortices[1].position[0] = DOMAIN_SIZE_X/2 + offset;
	vortices[1].position[1] = DOMAIN_SIZE_Y/2 + offset;
	vortices[1].intensity = vortMagnitude;
	vortices[1].initStep = 0;

	vortices[2].vIndex = 2;
	vortices[2].velocity = malloc(sizeof(double) * 2);
	vortices[2].position = malloc(sizeof(double) * 2);
	vortices[2].position[0] = DOMAIN_SIZE_X/2 - offset;
	vortices[2].position[1] = DOMAIN_SIZE_Y/2 - offset;
	vortices[2].intensity = vortMagnitude;
	vortices[2].initStep = 0;
	
	vortices[3].vIndex = 3;
	vortices[3].velocity = malloc(sizeof(double) * 2);
	vortices[3].position = malloc(sizeof(double) * 2);
	vortices[3].position[0] = DOMAIN_SIZE_X/2 + offset;
	vortices[3].position[1] = DOMAIN_SIZE_Y/2 - offset;
	vortices[3].intensity = vortMagnitude;
	vortices[3].initStep = 0;
}

void initialize_single_point(struct Vortex *vortices, int n) {
	vortices[0].vIndex = 0;
	vortices[0].velocity = malloc(sizeof(double) * 2);
	vortices[0].position = malloc(sizeof(double) * 2);
	vortices[0].position[0] = DOMAIN_SIZE_X/2;
	vortices[0].position[1] = DOMAIN_SIZE_Y/2;
	vortices[0].intensity = VORTEX_INTENSITY_INIT_UPPER_BOUND;
	vortices[0].initStep = 0;
}

void initialize_vorts(struct Vortex *vortices, int n) {
	switch (TEST_CASE) {
		case 0: {
			initialize_drivers_random(vortices, n, 0);
			break;
		} case 1: {
			initialize_pair_orbit_test(vortices, n);
			break;
		} case 2: {
			initialize_pair_parallel_test(vortices, n);
			break;
		} case 3: {
			initialize_square_system(vortices, n);
			break;
		} case 4: {
			initialize_single_point(vortices, n);
		}
		default:
			break;
	}
}

void initialize_tracers(struct Tracer *tracers, int n) {
	double seperationX = (DOMAIN_SIZE_X) / (sqrt(n) + 1);
	double seperationY = (DOMAIN_SIZE_Y) / (sqrt(n) + 1);

	int i = 0;

	for (int row = 1; row <= sqrt(n); row++) {
		for (int col = 1; col <= sqrt(n); col++) {
			tracers[i].tIndex = i;
			tracers[i].position = malloc(sizeof(double) * 2);
			tracers[i].velocity = calloc(sizeof(double), 2);

			tracers[i].position[0] = col * seperationX;
			tracers[i].position[1] = row * seperationY;
			i++;
		}
	}
}

#pragma mark - misc

void pprintVortRads(double *rads, int numActiveVorts) {
	long index = 0;
	for (int row = 0; row < numActiveVorts; row++) {
		for (int col = 0; col < row; col++) {
			printf("|");
			for (int i = 0; i < 3; i++) {
				printf("%6.2f", rads[index]);
				if (i != 2) printf(",");
				index++;
			}
			if (col == row-1) printf("|");
		}
		printf("\n");
	}
}

void pprintTracerRads(double *rads, int numActiveTracers, int numActiveVorts) {
	long index = 0;
	for (int vortIndex = 0; vortIndex < numActiveVorts; vortIndex++) {
		for (int tracerIndex = 0; tracerIndex < numActiveTracers; tracerIndex++) {
			printf("|");
			for (int i = 0; i < 3; i++) {
				printf("%6.2f", rads[index]);
				if (i != 2) printf(",");
				index++;
			}
		}
		printf("|\n");
	}
}

void wrapPositions(struct Vortex *vorts, int numVorts, struct Tracer *tracers, int numTracers) {
	for (int i = 0; i < numVorts; i++) {
		if (vorts[i].position[0] < 0) {
			vorts[i].position[0] = DOMAIN_SIZE_X + fmod(vorts[i].position[0], DOMAIN_SIZE_X);
		} else if (vorts[i].position[0] > DOMAIN_SIZE_X) {
			vorts[i].position[0] =  fmod(vorts[i].position[0], DOMAIN_SIZE_X);
		}
		if (vorts[i].position[1] < 0) {
			vorts[i].position[1] = DOMAIN_SIZE_Y + fmod(vorts[i].position[1], DOMAIN_SIZE_Y);
		} else if (vorts[i].position[1] > DOMAIN_SIZE_Y) {
			vorts[i].position[1] = fmod(vorts[i].position[1], DOMAIN_SIZE_Y);
		}
	}
	
	for (int i = 0; i < numTracers; i++) {
		if (tracers[i].position[0] < 0) {
			tracers[i].position[0] = DOMAIN_SIZE_X + fmod(tracers[i].position[0], DOMAIN_SIZE_X);;
		} else if (tracers[i].position[0] > DOMAIN_SIZE_X) {
			tracers[i].position[0] = fmod(tracers[i].position[0], DOMAIN_SIZE_X);
		}
		if (tracers[i].position[1] < 0) {
			tracers[i].position[1] = DOMAIN_SIZE_Y + fmod(tracers[i].position[1], DOMAIN_SIZE_Y);
		} else if (tracers[i].position[1] > DOMAIN_SIZE_Y) {
			tracers[i].position[1] = fmod(tracers[i].position[1], DOMAIN_SIZE_Y);
		}
	}
}

#pragma mark - main

float timespentDrawing = 0;

int main(int argc, const char * argv[]) {
	
	if (FIRST_SEED == -1) {
		time((time_t *)&randomSeed);
		printf("First random seed is %i\n", randomSeed);
	} else {
		randomSeed = FIRST_SEED;
	}

	currentTimestep = 0;
	int vorticesAllocated = (int)NUM_VORT_INIT*1.5;
	int activeDriverVortices = NUM_VORT_INIT;

	// vortices is the array of Vortex structs
	struct Vortex *vortices = malloc(sizeof(struct Vortex) * vorticesAllocated);
	initialize_vorts(vortices, activeDriverVortices);

	// radii is the matrix of distances between vortices. The distance between vortex vIndex==a and vortex vIndex==b (where a < b) is at index 3*(a*(a+1)/2+b).
	// the next item in the array is the x-component of the distance, and then the y-component of the distance
	// r, r_x, r_y

//	int radiiLen = sizeof(double) * (pow(vorticesAllocated, 2) - vorticesAllocated);
	unsigned long vortexRadiiLen = sizeof(double) * (calculateVortexRadiiIndex(vorticesAllocated-1, vorticesAllocated-2) + 1);

	double *vortexRadii = malloc(vortexRadiiLen);
//	double *radii = calloc(sizeof(double), pow(vorticesAllocated, 2) - vorticesAllocated);

	// calculate the radii for the first timestep. Update radii above only works after velocities are calculated, which we can't do until radii are calculated for the first time.

	long index;
	for (int i = 0; i < activeDriverVortices; ++i) {
		for (int j = 0; j < i; ++j) {
			index = calculateVortexRadiiIndex(vortices[i].vIndex, vortices[j].vIndex);
			vortexRadii[index+1] = vortices[i].position[0] - vortices[j].position[0];
			vortexRadii[index+2] = vortices[i].position[1] - vortices[j].position[1];
			vortexRadii[index] = sqrt(pow(vortexRadii[index+1], 2) + pow(vortexRadii[index+2], 2));

			double vTot = velocityFunc(vortices[j].intensity, vortexRadii[index]);
			vortices[i].velocity[0] = vTot*(vortexRadii[index+2]/vortexRadii[index]);
			vortices[i].velocity[1] = vTot*(-vortexRadii[index+1]/vortexRadii[index]);
		}
	}

	// tracers is the array of Tracer structs
	struct Tracer *tracers = malloc(sizeof(struct Tracer) * NUM_TRACERS);
	initialize_tracers(tracers, NUM_TRACERS);
	long tracerRadSize = NUM_TRACERS * vorticesAllocated * sizeof(double) * 3;
	double *tracerRadii = malloc(tracerRadSize); // row = vortex, col = tracer; Note: vortPos - tracerPos

	for (int tracerIndex = 0; tracerIndex < NUM_TRACERS; tracerIndex++) {
		for (int vortIndex = 0; vortIndex < activeDriverVortices; vortIndex++) {
			long radIndex = calculateTracerRadiiIndex(tracerIndex, vortIndex);
			tracerRadii[radIndex+1] = vortices[vortIndex].position[0] - tracers[tracerIndex].position[0];
			tracerRadii[radIndex+2] = vortices[vortIndex].position[1] - tracers[tracerIndex].position[1];
			tracerRadii[radIndex] = sqrt(pow(tracerRadii[radIndex+1], 2) + pow(tracerRadii[radIndex+2], 2));
		}
	}

//	pprintVortRads(vortexRadii, activeDriverVortices);
//	deleteVortex(vortices + 3, vortexRadii, &activeDriverVortices, vortices);
//	printf("----------------\n");
//	pprintVortRads(vortexRadii, activeDriverVortices);
	
	/******************* main loop *******************/

	while (NUMBER_OF_STEPS == 0 || currentTimestep < NUMBER_OF_STEPS) {
		struct timespec startTime;
		struct timespec endTime;
		
// #undef DEBUG
#ifdef DEBUG
//		printf("Current min r: %f\n", minR);
		if (NUMBER_OF_STEPS != 0 && !(currentTimestep%(NUMBER_OF_STEPS/20))) {
//		if (currentTimestep%RENDER_NTH_STEP == 0) {
			printf("-----------\n");
			for (int ind = 0; ind < NUM_TRACERS; ind += 1000) printf("%i	|	%1.5f	|	%1.5f	|	%1.5f\n",
																			  ind,
																			  tracers[ind].velocity[0],
																			  tracers[ind].velocity[1],
																			  sqrt(pow(tracers[ind].velocity[0], 2) + pow(tracers[ind].velocity[1], 2)));

			printf("Step number %i\n", currentTimestep);
			clock_gettime(CLOCK_MONOTONIC, &endTime);
			double sec = (endTime.tv_sec - startTime.tv_sec) + (double)(endTime.tv_nsec - startTime.tv_nsec) / 1E9;
			printf("Step took %f sec\n", sec);
		}
#else
#ifdef DRAW_CONSOLE
		if (currentTimestep%RENDER_NTH_STEP == 0) {
			drawToConsole(vortices, activeDriverVortices, tracers);
			for (int c = 0; c < ceil(log10(currentTimestep)) + 1; c++) printf("\010");
			printf("%i", currentTimestep);
			struct timespec sleepDuration;
			sleepDuration.tv_sec = 0;
			sleepDuration.tv_nsec = (1./5)*1E9;

			nanosleep(&sleepDuration, NULL);
		}
#endif
#ifdef DRAW_PDF
		if (currentTimestep%RENDER_NTH_STEP == 0) {
			char *filename = malloc(sizeof(char) * 50);
			genFName(filename, currentTimestep);
			
			clock_gettime(CLOCK_MONOTONIC, &startTime);
			drawToFile(vortices, activeDriverVortices, tracers, filename);
			clock_gettime(CLOCK_MONOTONIC, &endTime);
			timespentDrawing += (endTime.tv_sec - startTime.tv_sec) + (double)(endTime.tv_nsec - startTime.tv_nsec) / 1E9;
			
			
			printf("Finished frame: %s\n", filename);
		}
#endif
#endif
		
		stepForward_RK4(vortices, vortexRadii, activeDriverVortices, tracerRadii, tracers, NUM_TRACERS);
		wrapPositions(vortices, activeDriverVortices, tracers, NUM_TRACERS);
		updateRadii_pythagorean(vortexRadii, vortices, activeDriverVortices, tracerRadii, tracers, NUM_TRACERS); // needs optimization. Repeat ops from stepForward_RK4
		//		double minR = minRad(vortexRadii, activeDriverVortices);
		//		printf("Step number %i calculation complete\n", currentTimestep);
		currentTimestep++;
	}

	printf("Total time spent drawing: %5.2f sec\n", timespentDrawing);

	return 0;
}
