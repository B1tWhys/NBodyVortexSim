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

time_t test;

void testNormal() {
	int n = 100000;
	float resultsArr[n];
//	int totalSpawned = 0;
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
//	srand((unsigned int)time(0));
//	testNormal();
//	testPoisson();
	
	
	int var1 = 0;
	int var2 = 1;
	int var3 = 2;
	// bool __atomic_compare_exchange_n (type *ptr, type *expected, type desired, bool weak, int success_memorder, int failure_memorder)
	__atomic_compare_exchange_n(&var1, &var2, var3, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST);
	
	printf("1: %i\n2: %i\n3: %i\n", var1, var2, var3);
	
	return 0;
}
















