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

time_t test;


double generateUniformRandInRange(double lowerBound, double upperBound) { // range is inclusive on both ends
	// TODO: investigate whether there's an off by 1 error in here
	double range = upperBound-lowerBound;
	double result = ((double)rand()/RAND_MAX)*range + lowerBound;
	return result;
}


double nextNum = NAN;

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

void testFunc() {
	int n = 100000;
	float resultsArr[n];
	int totalSpawned = 0;
	for (int i = 0; i < n; i++) resultsArr[i] = 0;

	for (int i = 0; i < n; i++) {
		float val = generateNormalRand(1.);
		resultsArr[i] = val;
	}

	for (float cutoff = -5; cutoff < 5; cutoff++) {
		int count = 0;
		for (int i = 0; i < n; i++) {
			if (resultsArr[i] > cutoff && resultsArr[i] < (cutoff + 1)) count++;
		}
		printf("%i between %2.0f and %2.0f\n", count, cutoff, cutoff+1);
	}
}

#pragma mark - main

int main(int argc, const char * argv[]) {
	srand((unsigned int)time(0));
	testFunc();
	
	return 0;
}
















