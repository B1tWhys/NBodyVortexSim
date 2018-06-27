//
//  main.c
//  NBodySim
//
//  Created by Skyler Arnold on 6/12/18.
//  Copyright © 2018 Skyler Arnold. All rights reserved.
//

#include "main.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "constants.h"

#define BPPOINT printf("")

unsigned int randomSeed;
int currentTimestep;

long calculateVortexRadiiIndex(long vortIndex1, long vortIndex2) {
	// radii is basically a hash table, with this as the hash function.
	// I suspect that this is better than the standard dictionary for this.
	// I'm able to guarentee that the structure is full, and that there are no
	// collisions, and I'm able to iterate efficiently by incrementing an array pointer.

	if (vortIndex1 < vortIndex2) {
		return (vortIndex1*(vortIndex1+1) / 2) + vortIndex2 - 1;
	} else {
		return (vortIndex2*(vortIndex2+1) / 2) + vortIndex1 - 1;
	}
}

long calculateTracerRadiiIndex(long tracerIndex, long vortIndex) {
	return vortIndex * NUM_TRACERS + tracerIndex;
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

// This works reliably but it *may* be *very slightly* slower
void updateRadii_pythagorean(double *radii, struct Vortex *vortices, int numvortices) {
	long index;

	for (int i = 0; i < numvortices; ++i) { // if this doesn't segfault it'll be a goddamn miracle
		for (int j = 0; j < i; ++j) {
			index = calculateVortexRadiiIndex(vortices[i].vIndex, vortices[j].vIndex);

			radii[index+1] = vortices[i].position[0] - vortices[j].position[0];
			radii[index+2] = vortices[i].position[1] - vortices[j].position[1];
			radii[index] = sqrt(pow(radii[index+1], 2) + pow(radii[index+2], 2));
			++index;
		}
	}
}

double velocityFunc(double vortex2Intensity, double radius) { // calculates velocity for vortex 1 as a result of vortex 2
	return vortex2Intensity/(2*M_PI*radius);
}

void calculateDxDj_DyDj_vortex(double *dxdj, double *dydj, struct Vortex *vort, struct Vortex *vortices, int numVorts, double *intRads, long numIntRads) {
	for (int j = 0; j < numVorts; ++j) {
		if (vort->vIndex == j) {
			continue;
		}

		long radiiIndex = calculateVortexRadiiIndex(vort->vIndex, vortices[j].vIndex);
		long intRadiiIndex = radiiIndex;

		for (int domain = 0; domain <= 8; domain++) {
			double rad;
			double xRad = intRads[intRadiiIndex + 1], yRad = intRads[intRadiiIndex + 2];

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
						yRad += DOMAIN_SIZE_Y;
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

			if (rad > DOMAIN_SIZE_X) continue; // domain truncation

			double vmag = velocityFunc(vortices[j].intensity, rad); // total velocity magnitude
			*dxdj +=  (yRad/rad) * TIMESTEP * vmag;
			*dydj += (-xRad/rad) * TIMESTEP * vmag;
		}

		/*
		 We only store the x and y distances between each pair of vorts once per pair, so if this calculation is being done
		 backwards (we stored the radii for a->b but we're solving the influence of a->b) then the distance components are
		 the negative of what we want (we want b.x-a.x but what we've got is a.x-b.x). Whether its forward/backward depends on
		 the vIndex of each one.
		 */
		if (vort->vIndex > vortices[j].vIndex) {
			*dxdj = -*dxdj;
			*dydj = -*dydj;
}
	}
}

void calculateDxDj_DyDj_tracer(double *dxdj, double *dydj, struct Tracer *tracer, struct Tracer *tracers, int numTravers, double *intRads, long numIntRads, struct Vortex *vortices, int numVorts) {
	for (int vortIndex = 0; vortIndex < numVorts; vortIndex++) {
		long intRadIndex = calculateTracerRadiiIndex(tracer->tIndex, vortIndex);

		for (int domain = 0; domain <= 8; domain++) {
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
						yRad += DOMAIN_SIZE_Y;
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

			if (rad > DOMAIN_SIZE_X) continue;

			double vmag = velocityFunc(vortices[vortIndex].intensity, rad);
			*dxdj +=  (yRad/rad) * TIMESTEP * vmag;
			*dydj += (-xRad/rad) * TIMESTEP * vmag;
		}
	}
}

void stepForward_RK4(struct Vortex *vortices, double *vortRadii, int numVortices, double *tracerRadii, struct Tracer *tracers, int numTracers) {
	/*
	 its so inefficient     ୧(ಠ Д ಠ)୨

	 // this algorithm ends up making 2 coppies of each radius. One in each direction.
	 // TODO: I think I really want to just copy the whole damn radius matrix

	 // upon further consideration this was very stupid. Definately just do the whole matrix

	 *also im not sure if i am allowed to do this summing of k terms...?
	 */

	int sizeOfRadEntry = sizeof(double) * 3;
	long intermediateRadLength = (pow((long)numVortices, 2) - numVortices) * sizeOfRadEntry;
	double *intermediateVortRads = malloc(intermediateRadLength);
	memcpy(intermediateVortRads, vortRadii, intermediateRadLength);

	long intermediateTracerRadLength = numVortices * numTracers * sizeOfRadEntry;
	double *intermediateTracerRads = malloc(intermediateRadLength);
	memcpy(intermediateTracerRads, tracerRadii, intermediateTracerRadLength);

	long intermedateTracerRadLength = numTracers * numVortices * sizeOfRadEntry;

	for (int i = 0; i < numVortices; ++i) {
		struct Vortex *vort = &vortices[i];
		vort->velocity[0] = 0;
		vort->velocity[1] = 0;
	}


	for (int RKStep = 1; RKStep <= 4; RKStep++) { // NOTE: this adds a const (not an order) to the big-O of the alg
		double dxdj = 0, dydj = 0;

		for (int originVortIndex = 0; originVortIndex < numVortices; ++originVortIndex) {
			double k1_x = 0, k2_x = 0, k3_x = 0, k4_x = 0;
			double k1_y = 0, k2_y = 0, k3_y = 0, k4_y = 0;

			struct Vortex *vort = &vortices[originVortIndex];
			int vortIndex = vort->vIndex;

			int offset = 0; // neccessary for skipping over the self<->self radius in the vortices array

			dxdj = 0;
			dydj = 0;

			calculateDxDj_DyDj_vortex(&dxdj, &dydj, vort, vortices, numVortices, intermediateVortRads, numVortices);

			for (int j = 0; j < numVortices; ++j) {
				if (originVortIndex == j) {
					offset = -1;
					continue;
				}

				long radiiIndex = calculateVortexRadiiIndex(vortIndex, vortices[j].vIndex);
				long intRadiiIndex = radiiIndex;

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
			int tracerIndex = tracer->tIndex;

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
				long radIndex = (long)tracerIndex * (long)vortexIndex;

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

				intermediateVortRads[radIndex+1] = vortRadii[radIndex + 1] + dxdj;
				intermediateVortRads[radIndex+2] = vortRadii[radIndex + 2] + dydj;

				// pythagorean.
				intermediateVortRads[radIndex] = sqrt(pow(intermediateVortRads[radIndex+1], 2) + pow(intermediateVortRads[radIndex+2], 2));
			}
			tracer->velocity[0] += (k1_x + k2_x * 2 + k3_x * 2 + k4_x)/6;
			tracer->velocity[1] += (k1_y + k2_y * 2 + k3_y * 2 + k4_y)/6;
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

void loopVorts(struct Vortex *vorts, int n) {
	for (int i = 0; i < n; i++) {
		if (vorts[i].position[0] < 0) {
			vorts[i].position[0] = DOMAIN_SIZE_X;
		} else if (vorts[i].position[0] > DOMAIN_SIZE_X) {
			vorts[i].position[0] = 0;
		}
		if (vorts[i].position[1] < 0) {
			vorts[i].position[1] = DOMAIN_SIZE_Y;
		} else if (vorts[i].position[1] > DOMAIN_SIZE_Y) {
			vorts[i].position[1] = 0;
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
		vortices[i].intensity = generateRandInRange(0, VORTEX_INTENSITY_INIT_UPPER_BOUND);
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
	vortices[0].position[0] = DOMAIN_SIZE_X/2 - 10;
	vortices[0].position[1] = DOMAIN_SIZE_Y/2;
	vortices[0].intensity = intensity;
	vortices[0].initStep = 0;

	vortices[1].vIndex = 1;
	vortices[1].velocity = malloc(sizeof(double) * 2);
	vortices[1].position = malloc(sizeof(double) * 2);
	vortices[1].position[0] = DOMAIN_SIZE_X/2 + 10;
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
	vortices[0].position[0] = DOMAIN_SIZE_X/2 - 10;
	vortices[0].position[1] = DOMAIN_SIZE_Y/2;
	vortices[0].intensity = vortMagnitude;
	vortices[0].initStep = 0;

	vortices[1].vIndex = 1;
	vortices[1].velocity = malloc(sizeof(double) * 2);
	vortices[1].position = malloc(sizeof(double) * 2);
	vortices[1].position[0] = DOMAIN_SIZE_X/2 + 10;
	vortices[1].position[1] = DOMAIN_SIZE_Y/2;
	vortices[1].intensity = -vortMagnitude;
	vortices[1].initStep = 0;
}

void initialize_ternary_system(struct Vortex *vortices, int n) {
	assert(NUM_VORT_INIT == 3);
	assert(n == NUM_VORT_INIT);

	double vortMagnitude = VORTEX_INTENSITY_INIT_UPPER_BOUND;

	vortices[0].vIndex = 0;
	vortices[0].velocity = malloc(sizeof(double) * 2);
	vortices[0].position = malloc(sizeof(double) * 2);
	vortices[0].position[0] = DOMAIN_SIZE_X/2 - 10;
	vortices[0].position[1] = DOMAIN_SIZE_Y/2;
	vortices[0].intensity = vortMagnitude;
	vortices[0].initStep = 0;

	vortices[1].vIndex = 1;
	vortices[1].velocity = malloc(sizeof(double) * 2);
	vortices[1].position = malloc(sizeof(double) * 2);
	vortices[1].position[0] = DOMAIN_SIZE_X/2 + 10;
	vortices[1].position[1] = DOMAIN_SIZE_Y/2;
	vortices[1].intensity = vortMagnitude;
	vortices[1].initStep = 0;

	vortices[2].vIndex = 2;
	vortices[2].velocity = malloc(sizeof(double) * 2);
	vortices[2].position = malloc(sizeof(double) * 2);
	vortices[2].position[0] = DOMAIN_SIZE_X/2;
	vortices[2].position[1] = DOMAIN_SIZE_Y/2 + sqrt(3.0)*10;
	vortices[2].intensity = vortMagnitude;
	vortices[2].initStep = 0;
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
			initialize_ternary_system(vortices, n);
			break;
		}
		default:
			break;
	}
}

void initialize_tracers(struct Tracer *tracers, int n) {
	double seperationX = DOMAIN_SIZE_X / sqrt(n+2);
	double seperationY = DOMAIN_SIZE_Y / sqrt(n+2);

	for (int row = 1; row < sqrt(n) - 1; row++) {
		for (int col = 1; sqrt(col) < n - 1; col++) {
			int i = sqrt(n)*row + col;

			tracers[i].tIndex = i;
			tracers[i].position = malloc(sizeof(double) * 2);
			tracers[i].velocity = calloc(sizeof(double), 2);

			tracers[i].position[0] = (col + 1) * seperationX;
			tracers[i].position[1] = (row + 1) * seperationY;
		}
	}
}

#pragma mark - UI/Output code

void drawToConsole(struct Vortex *vorts, int n) {
	printf("\033[2J");
	int pxArray[CONSOLE_W][CONSOLE_H];

	for (int y = 0; y < CONSOLE_H; y++) {
		for (int x = 0; x < CONSOLE_W; x++) {
			pxArray[x][y] = -1;
		}
	}

	for (int i = 0; i < n; i++) {
		int xPxCoord = (int)(CONSOLE_W-1)*vorts[i].position[0]/DOMAIN_SIZE_X;
		int yPxCoord = (int)(CONSOLE_H-1)*vorts[i].position[1]/DOMAIN_SIZE_Y;

		pxArray[xPxCoord][yPxCoord] = vorts[i].vIndex;
	}
	int count = 0;
	for (int y = CONSOLE_H; y >= 0; y--) {
		for (int x = 0; x < CONSOLE_W; x++) {
			if (pxArray[x][y] != -1) {
				count++;
				for (int j = 0; j < ceil(log10(pxArray[x][y])); j++) printf("\010");
				printf("%i", pxArray[x][y]);
			} else {
				printf(" ");
			}
		}
		if (y != CONSOLE_H) printf("\n");
	}
}

#pragma mark - main

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

	// tracers is the arary of Tracer structs
	struct Tracer *tracers = malloc(sizeof(struct Tracer) * NUM_TRACERS);
	long tracerRadLength = NUM_TRACERS * vorticesAllocated * sizeof(double) * 3;
	double *tracerRadii = malloc(tracerRadLength); // row = vortex, col = tracer; Note: vortPos - tracerPos

	for (int tracerIndex = 0; tracerIndex < NUM_TRACERS; tracerIndex++) {
		for (int vortIndex = 0; vortIndex < activeDriverVortices; vortIndex++) {
			long radIndex = calculateTracerRadiiIndex(tracerIndex, vortIndex);
			tracerRadii[radIndex+1] = vortices[vortIndex].position[0] - tracers[tracerIndex].position[0];
			tracerRadii[radIndex+2] = vortices[vortIndex].position[1] - tracers[tracerIndex].position[1];
			tracerRadii[radIndex] = sqrt(pow(tracerRadii[radIndex+1], 2) + pow(tracerRadii[radIndex+2], 2));
		}
	}

	/******************* main loop *******************/

	while (NUMBER_OF_STEPS == 0 || currentTimestep < NUMBER_OF_STEPS) {
		struct timespec startTime;
		clock_gettime(CLOCK_MONOTONIC, &startTime);

		stepForward_RK4(vortices, vortexRadii, activeDriverVortices, tracerRadii, tracers, NUM_TRACERS);

		updateRadii_pythagorean(vortexRadii, vortices, activeDriverVortices);
		loopVorts(vortices, activeDriverVortices);

#ifdef DEBUG
		printf("Step number %i\n", currentTimestep);
		if (NUMBER_OF_STEPS != 0 && !(currentTimestep%(NUMBER_OF_STEPS/20))) {
			printf("-----------\n");
			for (int ind = 0; ind < activeDriverVortices; ind += 1000) printf("%i	|	%1.5f	|	%1.5f	|	%1.5f\n", ind, vortices[ind].velocity[0], vortices[ind].velocity[1], sqrt(pow(vortices[ind].velocity[0], 2) + pow(vortices[ind].velocity[1], 2)));

			printf("Step number %i\n", currentTimestep);
			struct timespec endTime;
			clock_gettime(CLOCK_MONOTONIC, &endTime);
			double sec = (endTime.tv_sec - startTime.tv_sec) + (double)(endTime.tv_nsec - startTime.tv_nsec) / 1E9;
			printf("Step took %f sec\n", sec);
		}
#else
		if (currentTimestep%RENDER_NTH_STEP == 0) {
			drawToConsole(vortices, activeDriverVortices);
			for (int c = 0; c < ceil(log10(currentTimestep)) + 1; c++) printf("\010");
			printf("%i", currentTimestep);
			struct timespec sleepDuration;
			sleepDuration.tv_sec = 0;
			sleepDuration.tv_nsec = .125E9;

			nanosleep(&sleepDuration, NULL);
		}
#endif
		currentTimestep++;
	}
}
