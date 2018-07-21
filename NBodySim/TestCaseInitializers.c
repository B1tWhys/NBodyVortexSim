//
//  TestCaseInitializers.c
//  NBodySim
//
//  Created by Skyler Arnold on 7/20/18.
//

#include "TestCaseInitializers.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>


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
void initialize_test(struct Vortex *vortices, int n) {
	switch (TEST_CASE) {
		case 1: {
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
