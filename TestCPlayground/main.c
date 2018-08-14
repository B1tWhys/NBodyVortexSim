//
//  main.c
//  TestCPlayground
//
//  Created by Skyler Arnold on 6/13/18.
//  Copyright © 2018 Skyler Arnold. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <signal.h>
#include "RNG.h"
#include "constants.h"

time_t test;
FILE *f;

void testNormal() {
	int n = 1E6;
	float resultsArr[n];
	
	for (int i = 0; i < n; i++) resultsArr[i] = 0;
	
	for (int i = 0; i < n; i++) {
		float val = generateNormalRand(VORTEX_INTENSITY_SIGMA);
		fprintf(f, "%.15f,", val);
		resultsArr[i] = val;
	}
	
	for (float cutoff = -5; cutoff < 5; cutoff++) {
		int count = 0;
		for (int i = 0; i < n; i++) {
			if (resultsArr[i] > cutoff && resultsArr[i] < (cutoff + 1)) count++;
		}
		printf("%f in %2.0f < x < %2.0f\n", (float)count/(float)n, cutoff, cutoff+1);
	}
}

void testPoisson() {
	int n = 10000;
	int data[100];
	for (int i = 0; i < 100; i++) data[i] = 0;
	
	for (int i = 0; i < n; i++) {
		int val = generatePoissonRand(0, 2.56, 0);
		data[val]++;
		printf("%i ", val);
	}
	printf("\n");
	float sum = 0;
	for (int i = 0; i < 15; i++) {
		printf("%i: %i\n", i, data[i]);
		sum += i*data[i];
	}
	
	printf("avg: %f", sum/10000);
}

#pragma mark - main

int main(int argc, const char * argv[]) {
	f = fopen("/Users/Skyler/Developer/randTestOutput", "w");
	srand(time(0));
//	FILE *f = fopen("/Users/Skyler/Developer/randTestOutput", "w");
//
//	for (int i = 0; i < 1E6; i++) {
//		fprintf(f, "%.15f,", generateNormalRand(1));
//	}
	printf("sigm: %f\n", VORTEX_INTENSITY_SIGMA);
	testNormal();
	
	return 0;
}
















