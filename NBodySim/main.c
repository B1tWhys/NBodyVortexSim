//
//  main.c
//  NBodySim
//
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
double timestep = TIMESTEP_CONST;
int currentTimestep;
#pragma mark - Index Calculators
/**
 This function return the index in the vortex radii array for the radius between a pair of vortices
 
 @param vortIndex1 the index of one of the vortices. Can be found by examining the vIndex value of the vortex in question
 @param vortIndex2 the index of the other vortex
 
 @return Index of the radius value
 */

long calculateVortexRadiiIndex(long vortIndex1, long vortIndex2) {
	// radii is basically a hash table, with this as the hash function.
	// I suspect that this is better than the standard dictionary libraries for this problem.
	// I'm able to guarentee that the structure is full, that there are no
	// collisions, and I'm able to iterate efficiently by incrementing an array pointer.
	if (vortIndex1 < vortIndex2) {
		return ((vortIndex2-1)*vortIndex2/2 + vortIndex1) * 3;
	} else {
		return ((vortIndex1-1)*vortIndex1/2 + vortIndex2) * 3;
	}
}

/**
 Calculates the index in the tracer radii array for the radius between a specific tracer and a specific vortex
 
 @param tracerIndex the index of the tracer. Can be found by exmaining the tIndex of the tracer in question
 @param vortIndex the index of the vortex. Can be found by examining the vIndex of the vortex in question
 
 @return The index in the tracer radii array of the radius
 */

long calculateTracerRadiiIndex(long tracerIndex, long vortIndex) {
	/*
	 rows = tracers
	 cols = vorts
	 */

	return (vortIndex * NUM_TRACERS + tracerIndex)*3;
}

#pragma mark - Diff Eq code

/**
 Recalculates all radii in a vortex radii array and a tracer radii array. Does so by calculating the pythagorean theorem for every radius. The radii arrays which are passed are modified.
 
 @param vortexRadii a pointer to the beginning of the array of doubles representing the radii of the vortices
 @param vortices a pointer to the beginning of the array of pointers to Vortex structs
 @param numVortices number of vortices in the vortexRadii array (which should be recalculated. This should not incldued space which was allocated but has not yet had a vortex assigned to it
 @param tracerRadii a pointer to the beginning of the array of doubles representing the radii between the tracers and vortices
 @param tracers a pointer to the array of Tracer structs
 @param numTracers the number of tracers which should be recalculated
*/
void updateRadii_pythagorean(double *vortexRadii, struct Vortex *vortices, int numVortices, double *tracerRadii, struct Tracer *tracers, int numTracers) {
	long index;

	for (int i = 0; i < numVortices; ++i) { // if this doesn't segfault it'll be a goddamn miracle
		for (int j = 0; j < i; ++j) {
			index = calculateVortexRadiiIndex(vortices[i].vIndex, vortices[j].vIndex);

			vortexRadii[index+1] = vortices[i].position[0] - vortices[j].position[0];
			vortexRadii[index+2] = vortices[i].position[1] - vortices[j].position[1];
			vortexRadii[index] = sqrt(pow(vortexRadii[index+1], 2.) + pow(vortexRadii[index+2], 2.));
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
			tracerRadii[index] = sqrt(pow(tracerRadii[index+1], 2.) + pow(tracerRadii[index+2], 2.));
		}
	}
}

/**
 Calculate the strength of the interaction of one body on a second body as a function of the distance between the current body and the 2nd one.
 
 @param vortex2Intensity The intensity of the other vortex
 @param radius The radial distance to the other vortex
 @return the magnitude of the velocity vector resulting from the interaction of the other vortex on the current vortex
 */
double velocityFunc(double vortex2Intensity, double radius) {
	return vortex2Intensity/(2.*M_PI*radius);
}

const char domains = 8; // 0 to disable wrapping of forces, 8 to enable it.
/**
 Calculate the x and y components of the velocity of a vortex based on the strengths of the other vortices, and the distance data in the intRads array
 @param dxdj Pointer to a double which will be updated to reflect the new x-velocity of the vortex
 @param dydj Pointer to a double which will be updated to reflect the new y-velocity of the vortex
 @param vort Pointer to the vortex to compute velocities for
 @param vortices Pointer to an array of all of the vortices in the system
 @param numVorts Number of vortices which are currently in the system
 @param rads The array containing all of the distance information between vortices as doubles
 @param numRads The length of the rads array
 */
void calculateDxDj_DyDj_vortex(double *dxdj, double *dydj, struct Vortex *vort, struct Vortex *vortices, int numVorts, double *rads, long numRads) {
	for (int j = 0; j < numVorts; ++j) {
		if (vort->vIndex == j) {
			continue;
		}
		
		long radiiIndex = calculateVortexRadiiIndex(vort->vIndex, vortices[j].vIndex);

		for (int domain = 0; domain <= domains; domain++) {
			double rad;
			double xRad = (vort->vIndex < vortices[j].vIndex) ? rads[radiiIndex + 1] : -rads[radiiIndex + 1];
			double yRad = (vort->vIndex < vortices[j].vIndex) ? rads[radiiIndex + 2] : -rads[radiiIndex + 2];
			
			if (domain) {
				/* domain numbering scheme:
				 1 2 3
				 4 0 5
				 6 7 8
				 */
				
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
				rad = rads[radiiIndex];
			}
			
			if (rad > DOMAIN_SIZE_X) {
				continue; // domain truncation
			}
			
			double vmag = velocityFunc(vortices[j].intensity, rad); // total velocity magnitude
			*dxdj +=  (yRad/rad) * vmag;
			*dydj += (-xRad/rad) * vmag;
		}
	}
}

/**
 Calculate the x and y components of the velocity of a tracer based on the strengths of the vortices and their distances from the tracer
 @param dxdj Pointer to a double which will be updated to reflect the new x-velocity of the vortex
 @param dydj Pointer to a double which will be updated to reflect the new y-velocity of the vortex
 @param tracer Pointer to the tracer to compute velocities for
 @param numVorts Number of vortices which are currently in the system
 @param rads The array containing all of the distance information between tracers and vortices as doubles
 @param numRads The length of the rads array
 */
void calculateDxDj_DyDj_tracer(double *dxdj, double *dydj, struct Tracer *tracer, double *rads, long numRads, struct Vortex *vortices, int numVorts) {
	for (int vortIndex = 0; vortIndex < numVorts; vortIndex++) {
		long intRadIndex = calculateTracerRadiiIndex(tracer->tIndex, vortIndex);

		for (int domain = 0; domain <= domains; domain++) {
			double rad;
			double xRad = rads[intRadIndex + 1], yRad = rads[intRadIndex + 2];
			
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
				rad = rads[intRadIndex];
			}
			
			if (rad > DOMAIN_SIZE_X || rad < .1) {
				continue;
			}

			double vmag = velocityFunc(vortices[vortIndex].intensity, rad);
			*dxdj +=  (yRad/rad) * vmag;
			*dydj += (-xRad/rad) * vmag;
		}
	}
}

/**
 moves the simulation forward 1 timestep using runge-kutta 4th order. This function updates all vortex and tracer positions and velocities
 @note This function does not update the appropriate radius arrays. Call @c UpdateRadii_Pythagorean() to update those arrays.
 
 @param vortices The array of all of the vortices in the simulation
 @param vortRadii The array of doubles containing vortex <-> vortex distance information
 @param numVortices The number of vortices currently in the simulation
 @param tracerRadii The array of doubles containing vortex <-> tracer distance information
 @param numTracers The numebr of tracers
 */
void stepForward_RK4(struct Vortex *vortices, double *vortRadii, int numVortices, double *tracerRadii, struct Tracer *tracers, int numTracers) {
	/*
	 There may be some duplicate computation being done here.     ୧(ಠ Д ಠ)୨
	 */

	int sizeOfRadEntry = sizeof(double) * 3;
	long vortRadSize = (pow(numVortices, 2)-numVortices)/2 * sizeOfRadEntry;
	
	/*
	 workingRadii is updated after every vortex is updated in position, and should not be directly used to compute vortex velocities
	 intermediateRadii is updated after every runge-kutta step. At the end of the step, workingRadii is coppied into intermediateRadii.
	 	IntermediateRadii should be used to compute the velocities of each particle.
	 */
	double *workingRadii = malloc(vortRadSize);
	memcpy(workingRadii, vortRadii, vortRadSize);
	double *intermediateRadii = malloc(vortRadSize);
	memcpy(intermediateRadii, vortRadii, vortRadSize);
	
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
		
		for (int tracerIndex = 0; tracerIndex < numTracers; ++tracerIndex) {
			double k1_x = 0, k2_x = 0, k3_x = 0, k4_x = 0;
			double k1_y = 0, k2_y = 0, k3_y = 0, k4_y = 0;
			
			struct Tracer *tracer = &tracers[tracerIndex];
			dxdj = 0;
			dydj = 0;
			
			calculateDxDj_DyDj_tracer(&dxdj,
									  &dydj,
									  tracer,
									  intermediateTracerRads,
									  (numVortices * numTracers),
									  vortices,
									  numVortices);
			
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
			
#ifdef DEBUG
			switch (RKStep) {
				case 1: {
					printf("step: %i | tracer: %i | k1_x: %1.15f\n", currentTimestep, tracer->tIndex, k1_x);
					printf("step: %i | tracer: %i | k1_y: %1.15f\n", currentTimestep, tracer->tIndex, k1_y);
					break;
				}
				case 2: {
					printf("step: %i | tracer: %i | k2_x: %1.15f\n", currentTimestep, tracer->tIndex, k2_x);
					printf("step: %i | tracer: %i | k2_y: %1.15f\n", currentTimestep, tracer->tIndex, k2_y);
					break;
				}
				case 3: {
					printf("step: %i | tracer: %i | k3_x: %1.15f\n", currentTimestep, tracer->tIndex, k3_x);
					printf("step: %i | tracer: %i | k3_y: %1.15f\n", currentTimestep, tracer->tIndex, k3_y);
					break;
				}
				case 4: {
					printf("step: %i | tracer: %i | k4_x: %1.15f\n", currentTimestep, tracer->tIndex, k4_x);
					printf("step: %i | tracer: %i | k4_y: %1.15f\n", currentTimestep, tracer->tIndex, k4_y);
					break;
				}
				default: {
					break;
				}
			}
#endif
		for (int vortexIndex = 0; vortexIndex < numVortices; ++vortexIndex) {
			long radIndex = calculateTracerRadiiIndex(tracerIndex, vortexIndex);
			intermediateTracerRads[radIndex + 1] = tracerRadii[radIndex + 1] + dxdj * timestep;
			intermediateTracerRads[radIndex + 2] = tracerRadii[radIndex + 2] + dydj * timestep;
			
			// pythagorean.
			intermediateTracerRads[radIndex] = sqrt(pow(intermediateTracerRads[radIndex+1], 2) + pow(intermediateTracerRads[radIndex+2], 2));
		}
			tracer->velocity[0] += (k1_x + k2_x * 2. + k3_x * 2. + k4_x)/6.;
			tracer->velocity[1] += (k1_y + k2_y * 2. + k3_y * 2. + k4_y)/6.;
		}
		
		for (int originVortIndex = 0; originVortIndex < numVortices; ++originVortIndex) {
			double k1_x = 0, k2_x = 0, k3_x = 0, k4_x = 0;
			double k1_y = 0, k2_y = 0, k3_y = 0, k4_y = 0;
			
			struct Vortex *vort = &vortices[originVortIndex];
			int vortIndex = vort->vIndex;
			
			dxdj = 0;
			dydj = 0;
			
			calculateDxDj_DyDj_vortex(&dxdj, &dydj, vort, vortices, numVortices, intermediateRadii, numVortices);
			
			switch (RKStep) {
				case 1: {
					k1_x = dxdj;
					k1_y = dydj;
					
					dxdj = dxdj/2.;
					dydj = dydj/2.;
					break;
				}
				case 2: {
					k2_x = dxdj;
					k2_y = dydj;
					
					dxdj = dxdj/2.;
					dydj = dydj/2.;
					break;
				}
				case 3: {
					k3_x = dxdj;
					k3_y = dydj;
					break;
				}
				case 4: {
					k4_x = dxdj;
					k4_y = dydj;
					break;
				}
				default: {
					break;
				}
			}
			
#ifdef DEBUG
			switch (RKStep) {
				case 1: {
					printf("step: %i | vortex: %i | k1_x: %1.15f\n", currentTimestep, vort->vIndex, k1_x);
					printf("step: %i | vortex: %i | k1_y: %1.15f\n", currentTimestep, vort->vIndex, k1_y);
					break;
				}
				case 2: {
					printf("step: %i | vortex: %i | k2_x: %1.15f\n", currentTimestep, vort->vIndex, k2_x);
					printf("step: %i | vortex: %i | k2_y: %1.15f\n", currentTimestep, vort->vIndex, k2_y);
					break;
				}
				case 3: {
					printf("step: %i | vortex: %i | k3_x: %1.15f\n", currentTimestep, vort->vIndex, k3_x);
					printf("step: %i | vortex: %i | k3_y: %1.15f\n", currentTimestep, vort->vIndex, k3_y);
					break;
				}
				case 4: {
					printf("step: %i | vortex: %i | k4_x: %1.15f\n", currentTimestep, vort->vIndex, k4_x);
					printf("step: %i | vortex: %i | k4_y: %1.15f\n", currentTimestep, vort->vIndex, k4_y);
					break;
				}
				default: {
					break;
				}
			}
#endif
			
			for (int j = 0; j < numVortices; ++j) {
				if (originVortIndex == j) {
					continue;
				}
				
				long radiiIndex = calculateVortexRadiiIndex(vortIndex, vortices[j].vIndex);
				// there is probably a cleverer way to do this to avoid having to copy radii into workingRadii every timestep
				// Whether we add or subtract d_dj from radii depends on which index is larger.
				if (vort->vIndex < vortices[j].vIndex) {
					workingRadii[radiiIndex + 1] -= dxdj*timestep;
					workingRadii[radiiIndex + 2] -= dydj*timestep;
				} else {
					workingRadii[radiiIndex + 1] += dxdj*timestep;
					workingRadii[radiiIndex + 2] += dydj*timestep;
				}
				
				workingRadii[radiiIndex] = sqrt(pow(workingRadii[radiiIndex+1], 2) + pow(workingRadii[radiiIndex+2], 2));
			}
			
			for (int tracerI = 0; tracerI < numTracers; tracerI++) {
				long index = calculateTracerRadiiIndex(tracerI, originVortIndex);
				intermediateTracerRads[index + 1] -= dxdj*timestep;
				intermediateTracerRads[index + 2] -= dydj*timestep;
				intermediateTracerRads[index] = sqrt(pow(intermediateTracerRads[index + 1], 2) + pow(intermediateTracerRads[index + 2], 2));
			}
			
			vort->velocity[0] += (k1_x + k2_x * 2 + k3_x * 2 + k4_x)/6;
			vort->velocity[1] += (k1_y + k2_y * 2 + k3_y * 2 + k4_y)/6;
		}
		
		memcpy(intermediateRadii, workingRadii, vortRadSize);
		memcpy(workingRadii, vortRadii, vortRadSize);
	}
	
	for (int i = 0; i < numVortices; ++i) {
		struct Vortex *vort = &vortices[i];
		vort->position[0] += vort->velocity[0] * timestep;
		vort->position[1] += vort->velocity[1] * timestep;
	}
	
	for (int i = 0; i < numTracers; ++i) {
		struct Tracer *tracer = &tracers[i];
		tracer->position[0] += tracer->velocity[0] * timestep;
		tracer->position[1] += tracer->velocity[1] * timestep;
	}
	
	free(workingRadii);
	free(intermediateRadii);
	free(intermediateTracerRads);
}

/**
 generates a uniformly random double value in a range
 
 @param lowerBound the lower bound for the random number
 @param upperBound the upper bound for the random number
 
 @return returns a random uniform double in the given range
 */
double generateUniformRandInRange(double lowerBound, double upperBound) { // range is inclusive on both ends
	// TODO: investigate whether there's an off by 1 error in here
	double range = upperBound-lowerBound;
	double result = ((double)random()/RAND_MAX)*range + lowerBound;
	return result;
}

double nextNum = NAN;

/**
 generates a random double from a normal distributuion centered at 0 with a given standard deviation
 
 @discussion this uses the Box-Muller transform to create pairs of random numbers.
 	one is returned, and the other is stored for the next time that the function is called.
 
 @param sigma the standard deviation for the normal distribution
 */
double generateNormalRand(double sigma) {
	if (nextNum == NAN) {// if there is a number left over from last time, cache it, clear nextNum, and return the value.
		double val = nextNum;
		nextNum = NAN;
		return val;
	}
	
	double s, R, u, v;
	
	while (1) { // generate u & v until s is in (0, 1]
		u = generateUniformRandInRange(-1, 1);
		v = generateUniformRandInRange(-1, 1);
		
		s = pow(u, 2) + pow(v, 2);
		if (s != 0 && s < 1) break;
	}
	
	R = pow(s, 2);
	double multiplier = sqrt(-2.*log(s) / s);
	nextNum = u * multiplier;
	return v * multiplier;
}

#pragma mark - Designated Initializers

/**
 find the smallest radius seperation between vortices
 
 @param radArr the array of doubles contianing all vortex radius information
 @param numVorts the number of active driver vortices
 */
double minRad(double *radArr, long numVorts) {
	if (numVorts == 0) return 0;
	long radLen = (pow(numVorts, 2)-numVorts)/2 * 3;

	double min = radArr[0];
	for (int i = 0; i < radLen; i += 3) {
		if (radArr[i] < min) min = radArr[i];
	}
	return min;
}

/**
 delete a vortex and all associated data from the simulation
 
 @param vort the vortex to remove
 @param vortexRadii the array of doubles containing all vortex radius data
 @param numActiveVorts a pointer to the number of active vortices. The value pointed to will be decremented
 @param vorts the array of vortices
 @param tracerRads the array of doubles containing all tracer <-> vortex radius data
 */
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
/**
 find all vortices which should be merged and perform the merges
 
 @param vortexRadii the array of doubles containing all vortex radii data
 @param numActiveVorts a pointer to the number of active driver vorts. This value will be updated by the function to reflect merges
 @param vorts the array of vortices
 @param tracerRads the array of doubles containing all tracer <-> vortex radius data
 */
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

/**
 initialize the driver vortices (in the 1st domain) with and random strengths and positions
 
 @param vortices an empty array to put the new vortices in
 @param n number of vortices to initialize
 @param startingID the vID of the first new vortex
 */
void initialize_drivers_random(struct Vortex *vortices, int n, int startingID) {
	for (int i = 0; i < n; ++i) {
		vortices[i].vIndex = startingID++;
		vortices[i].position = malloc(sizeof(double) * 2);
		vortices[i].position[0] = generateUniformRandInRange(0, DOMAIN_SIZE_X);
		vortices[i].position[1] = generateUniformRandInRange(0, DOMAIN_SIZE_Y);
		vortices[i].velocity = calloc(2, sizeof(double));
		vortices[i].intensity = generateNormalRand(VORTEX_INTENSITY_SIGMA);
		vortices[i].initStep = 0;
	}
}

void initialize_pair_orbit_test(struct Vortex *vortices, int n) {
	assert(NUM_VORT_INIT == 2);
	assert(n == NUM_VORT_INIT);
	double seperation = 4;
	int intensity = 10;

	vortices[0].vIndex = 0;
	vortices[0].velocity = calloc(sizeof(double), 2);
	vortices[0].position = malloc(sizeof(double) * 2);
	vortices[0].position[0] = DOMAIN_SIZE_X/2 - seperation/2.;
	vortices[0].position[1] = DOMAIN_SIZE_Y/2;
	vortices[0].intensity = intensity;
	vortices[0].initStep = 0;

	vortices[1].vIndex = 1;
	vortices[1].velocity = calloc(sizeof(double), 2);
	vortices[1].position = malloc(sizeof(double) * 2);
	vortices[1].position[0] = DOMAIN_SIZE_X/2 + seperation/2.;
	vortices[1].position[1] = DOMAIN_SIZE_Y/2;
	vortices[1].intensity = intensity;
	vortices[1].initStep = 0;
}

void initialize_pair_parallel_test(struct Vortex *vortices, int n) {
	assert(NUM_VORT_INIT == 2);
	assert(n == NUM_VORT_INIT);

	double vortMagnitude = 5;

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

	double vortMagnitude = 1;
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
	vortices[0].intensity = 1;
	vortices[0].initStep = 0;
}

void initialize_test_case_4(struct Vortex *vortices, int n) {
	assert(n == 3);
	
	double l12 = 10;
	double thk = M_PI_4;
//	double tc = (5.0 - 3.0 * cos(2.0 * thk))/12.0/sin(2.0 * thk) * pow(l12, 2);
//	printf("Collapse time: %5.2f\n", tc);
	
	vortices[0].vIndex = 0;
	vortices[0].velocity = calloc(sizeof(double), 2);
	vortices[0].position = malloc(sizeof(double) * 2.);
	vortices[0].position[0] = (3. + sqrt(3.)*cos(thk)) / 6.0 * l12 + DOMAIN_SIZE_X / 2.;
	vortices[0].position[1] = sqrt(3.) * sin(thk) / 6.0 * l12 + DOMAIN_SIZE_Y / 2.;
	vortices[0].intensity = 4.*M_PI;
	vortices[0].initStep = 0;
	
	vortices[1].vIndex = 1;
	vortices[1].velocity = calloc(sizeof(double), 2);
	vortices[1].position = malloc(sizeof(double) * 2.);
	vortices[1].position[0] = (-3. + sqrt(3.)*cos(thk)) / 6.0 * l12 + DOMAIN_SIZE_X / 2.;
	vortices[1].position[1] = sqrt(3.) * sin(thk) / 6.0 * l12 + DOMAIN_SIZE_Y / 2.;
	vortices[1].intensity = 4.*M_PI;
	vortices[1].initStep = 0;
	
	vortices[2].vIndex = 2;
	vortices[2].velocity = calloc(sizeof(double), 2);
	vortices[2].position = malloc(sizeof(double) * 2.);
	vortices[2].position[0] = 2. * sqrt(3.) * cos(thk) / 3.0 * l12 + DOMAIN_SIZE_X / 2.;
	vortices[2].position[1] = 2. * sqrt(3.) * sin(thk) / 3.0 * l12 + DOMAIN_SIZE_Y / 2.;
	vortices[2].intensity = -2.*M_PI;
	vortices[2].initStep = 0;
}

/**
 handles initializeing vortices using the appropriate initializer for the current test case
 @param vortices array to put the new vortices into
 @param n the number of vortices to create
 */
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
			initialize_test_case_4(vortices, n);
			break;
		} case 5: {
			initialize_single_point(vortices, n);
			break;
		} case 6: {
			initialize_test_case_4(vortices, n);
			break;
		}
		default:
			break;
	}
}

/**
 generate the tracers
 
 @discussion this function does not put any tracers on the edge of the domain. If the space between tracers at t=0 is 1, there is a gap of size 1
 between each edge of the center domain and the closest column/row of tracers.
 
 @param tracers the array to fill with tracers
 @param n the number of tracers to generate
 */
void initialize_tracers(struct Tracer *tracers, int n) {
	double seperationX = (DOMAIN_SIZE_X) / (sqrt(n) + 1);
	double seperationY = (DOMAIN_SIZE_Y) / (sqrt(n) + 1);
	
	assert(sqrt(n) == (int)sqrt(n));
	
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

void initialize_single_test_tracer(struct Tracer *tracers, int numTracers, struct Vortex *vorts) {
	assert(numTracers == 1);
	
	initialize_tracers(tracers, numTracers);
	memcpy(tracers[0].position, vorts[0].position, sizeof(double) * 2);
}

#pragma mark - misc
/**print the vortRads array in a human readable format*/
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
/** print the tracerRads array in a human readable format*/
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
/**find vortices/tracers which have moved beyond the edges of the driver domain, and move them to the opposite side of the array.*/
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
/**
 find the highest velocity vortex
 
 @return the magnitude of the largest velocity vector
 */
double maxVelocity(struct Vortex *vorts, int numVorts) {
	double maxV = 0;
	for (int i = 0; i < numVorts; i++) {
		struct Vortex *vort = &vorts[i];
		double totalV = sqrt(pow(vort->velocity[0], 2) + pow(vort->velocity[1], 2));
		if (totalV > maxV) maxV = totalV;
	}
	return maxV;
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

//	int radiiLen = sizeof(double) * (pow(vorticesAllocated, 2) - vorticesAllocated);
	unsigned long vortexRadiiLen = sizeof(double) * (calculateVortexRadiiIndex(vorticesAllocated-1, vorticesAllocated-2) + 3);

	// vortexRadii is the matrix of distances between vortices. The distance between vortex vIndex==a and vortex vIndex==b (where a < b) is at index 3*(a*(a+1)/2+b).
	// the next item in the array is the x-component of the distance, and then the y-component of the distance
	// r, r_x, r_y
	
	double *vortexRadii = malloc(vortexRadiiLen);
//	double *radii = calloc(sizeof(double), pow(vorticesAllocated, 2) - vorticesAllocated);

	// tracers is the array of Tracer structs
	struct Tracer *tracers = malloc(sizeof(struct Tracer) * NUM_TRACERS);
	
	if (TEST_CASE != 6) {
		initialize_tracers(tracers, NUM_TRACERS);
	} else {
		initialize_single_test_tracer(tracers, NUM_TRACERS, vortices);
	}
	
	
	long tracerRadSize = NUM_TRACERS * vorticesAllocated * sizeof(double) * 3;
	double *tracerRadii = malloc(tracerRadSize); // row = vortex, col = tracer; Note: vortPos - tracerPos
	updateRadii_pythagorean(vortexRadii, vortices, activeDriverVortices, tracerRadii, tracers, NUM_TRACERS);

	/******************* main loop *******************/

	double time = 0;
	
	while (NUMBER_OF_STEPS == 0 || currentTimestep < NUMBER_OF_STEPS) {
		struct timespec startTime;
		struct timespec endTime;

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
		}
#endif
		// adjust the timestep based on minRad and maxVel for test case 4
		if (TEST_CASE == 4) {
			double minR = minRad(vortexRadii, activeDriverVortices);
			double maxV = maxVelocity(vortices, activeDriverVortices);
			timestep = minR / maxV * .5;
			if (timestep > TIMESTEP_CONST || maxV == 0) timestep = TIMESTEP_CONST;
			time += timestep;
			
			if (time > 50) return 0;
		}
		clock_gettime(CLOCK_MONOTONIC, &startTime);
		
		stepForward_RK4(vortices, vortexRadii, activeDriverVortices, tracerRadii, tracers, NUM_TRACERS);
		wrapPositions(vortices, activeDriverVortices, tracers, NUM_TRACERS);
		updateRadii_pythagorean(vortexRadii, vortices, activeDriverVortices, tracerRadii, tracers, NUM_TRACERS); // needs optimization. Repeat ops from stepForward_RK4
		
		clock_gettime(CLOCK_MONOTONIC, &endTime);
		double sec = (endTime.tv_sec - startTime.tv_sec) + (double)(endTime.tv_nsec - startTime.tv_nsec) / 1E9;
		printf("Step number %i calculation complete in %f sec\n", currentTimestep, sec);
		currentTimestep++;
	}

	printf("Total time spent drawing: %5.2f sec\n", timespentDrawing);

	return 0;
}

