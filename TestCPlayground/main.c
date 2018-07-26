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

int generatePoissonRand(double k, double lambda, double x) {
	double p, s, u;
	x = 0;
	p = pow(M_E, -lambda);
	s = p;
	
	u = generateUniformRandInRange(0, 1);
	
	while (u > s) {
		x++;
		p = (p * lambda) / x;
		s += p;
	}
	return x;
}

void testFunc() {
	int n = 100000;
	int distArr[100];
	int totalSpawned = 0;
	for (int i = 0; i < 100; i++) distArr[i] = 0;

	for (int i = 0; i < n; i++) {
		int val = generatePoissonRand(0, 1., 0);
		totalSpawned += val;
		distArr[val]++;
	}

	for (int i = 0; i < 10; i++) {
		printf("%i: %f\n", i, distArr[i]/(float)n);
	}

	printf("\n\ntotalSpawned: %i\n", totalSpawned);
}

//void testFunc() {
//	for (int i = 0; i < 5; i++) {
////		printf("%i\n", generatePoissonRand(0, 1., 0));
//		printf("%f\n", generateUniformRandInRange(0, 1));
//	}
//}

#pragma mark - main

int main(int argc, const char * argv[]) {
	test = time(0);
	srand((unsigned int)test);
	testFunc();
	
	return 0;
}
















