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

#pragma mark - test1
#define PRINT_DEBUG_N 1



int f(int x);
int g(int x);
int h(int x, int y);
void printFirstN(int *arr, int n);

void test1() {
	unsigned int ptrSize = sizeof(void *);
	unsigned int longSize = sizeof(void *);
	
	int numVals = 1000000;
	int n = numVals*3;
	
	int *dynamicArray = (int *)calloc(n, sizeof(int));
	
	struct timespec startTime, endTime;
	clock_gettime(CLOCK_MONOTONIC, &startTime);
	
	int i = 0;
	int j = 0;
	while (i < n) {
		int val1 = f(j);
		int val2 = g(j);
		++j;
		dynamicArray[i] = h(val1, val2);
		++i;
		dynamicArray[i] = val1;
		++i;
		dynamicArray[i] = val2;
		++i;
	}
	
	clock_gettime(CLOCK_MONOTONIC, &endTime);
	long ns = (endTime.tv_sec - startTime.tv_sec) * 1E9 + (endTime.tv_nsec - startTime.tv_nsec);
	printf("Test 1 took %li ns\n", ns);
	if (PRINT_DEBUG_N) printFirstN(dynamicArray, PRINT_DEBUG_N);
	
	
	free(dynamicArray);
	dynamicArray = (int *)calloc(n, sizeof(int));
	
	clock_gettime(CLOCK_MONOTONIC, &startTime);
	
	for (i = 0; i < numVals; ++i) {
		dynamicArray[i*3+1] = f(i);
		dynamicArray[i*3+2] = g(i);
		dynamicArray[i*3] = h(dynamicArray[i*3+1], dynamicArray[i*3+2]);
	}
	
	clock_gettime(CLOCK_MONOTONIC, &endTime);
	ns = (endTime.tv_sec - startTime.tv_sec) * 1E9 + (endTime.tv_nsec - startTime.tv_nsec);
	printf("Test 2 took %li ns\n", ns);
	if (PRINT_DEBUG_N) printFirstN(dynamicArray, PRINT_DEBUG_N);
	
	
	
	free(dynamicArray);
	dynamicArray = (int *)calloc(n, sizeof(int));
	
	clock_gettime(CLOCK_MONOTONIC, &startTime);
	
	i = 0;
	j = 0;
	while (i < n) {
		++i;
		dynamicArray[i] = f(j);
		++i;
		dynamicArray[i] = g(j);
		dynamicArray[i-2] = h(dynamicArray[i], dynamicArray[i-1]);
		++i;
		++j;
	}
	
	clock_gettime(CLOCK_MONOTONIC, &endTime);
	ns = (endTime.tv_sec - startTime.tv_sec) * 1E9 + (endTime.tv_nsec - startTime.tv_nsec);
	printf("Test 3 took %li ns\n", ns);
	if (PRINT_DEBUG_N) printFirstN(dynamicArray, PRINT_DEBUG_N);
}

void printFirstN(int *arr, int n) {
	printf("----------\n");
	for (int i = 0; i < n; i++) {
		printf("%2i: %i\n", i, arr[i]);
	}
	printf("----------\n");
}

int f(int x) {
	return x*2;
}

int g(int x) {
	return f(x+1);
}

int h(int x, int y) {
	return f(x+y);
}

#pragma mark - test 2

int randomSeed;

double generateRandInRange(double lowerBound, double upperBound) { // range is inclusive on both ends. if you flip the order you get a negative random #
	srand(randomSeed++);
	double range = upperBound-lowerBound;
	double result = (double)random()/(RAND_MAX)*(range) + lowerBound;
	return result;
}

void test2() {
	time((time_t *)&randomSeed);
	printf("seed 0: %i\n", randomSeed);
	
	for (int i = 0; i < 2000; i++) {
		double rand = generateRandInRange(0, 10);
		printf("%lf\n", rand);
	}
	
}

#pragma mark - test 3

struct testStruct1 {
	unsigned int *testArr;
};

void test3() {
	struct testStruct1 *testS = malloc(sizeof(struct testStruct1));
	testS->testArr = (int *)malloc(sizeof(int)*2);
	testS->testArr[0] = 1;
	testS->testArr[1] = 2;
	
	unsigned long val1 = 4*((unsigned long)UINT_MAX+1);
	unsigned long val2 = 3;
	
	unsigned long *doubleVal = malloc(sizeof(long));
	*doubleVal = val1 + val2;
	
	testS->testArr = (unsigned int *)doubleVal;
	printf("val1: %li\nval2: %li\n", testS->testArr[0], testS->testArr[1]);
}

#pragma mark - test 4

void test4() {
#ifdef TESTVAL
	printf("a: %i\n", TESTVAL);
#else
	printf("b\n");
#endif
}

#pragma mark - test 5

#define START_TIMER \
clock_gettime(CLOCK_MONOTONIC, &startTime);

#define STOP_TIMER \
clock_gettime(CLOCK_MONOTONIC, &endTime); \
sec = (endTime.tv_sec - startTime.tv_sec) + (double)(endTime.tv_nsec - startTime.tv_nsec) / 1E9;

void testFunc5() {
	//	unsigned char testBitwiseVal = 0;
	//
	//	testBitwiseVal += 0b00010000;
	//	testBitwiseVal += 0b01000000;
	//	testBitwiseVal += 0b00000100;
	//
	//	for (int i = 0; i < 8; i++) {
	//		if (testBitwiseVal & (unsigned char)pow(2, i)) {
	//			printf("%i bit had a value\n", i);
	//		}
	//	}
	
	struct timespec startTime, endTime;
	float sec;
	
	int len = 10000000;
	
	double *testArr1 = malloc(sizeof(double) * len);
	double *testArr2 = malloc(sizeof(double) * len);
	double *testArr3 = calloc(sizeof(double), len);
	
	
	for (int i = 0; i < len; i++) {
		testArr1[i] = (double)rand()/RAND_MAX;
		testArr2[i] = (double)rand()/RAND_MAX;
	}
	
	START_TIMER
	//	struct timespec startTime;
	//	clock_gettime(CLOCK_MONOTONIC, &startTime);
	
	for (int i = 0; i < len; i++) {
		testArr3[i] = testArr1[i] + testArr2[i];
	}
	
	STOP_TIMER
	printf("addition test took: %f sec\n", sec);
	printf("%f\n", testArr3[100]); // we have to print/use the contents of the array so that the compiler doesn't just remove the loop when optimiztaion ≠ lvl 0
	
	free(testArr3);
	testArr3 = calloc(sizeof(double), len);
	
	START_TIMER
	
	for (int i = 0; i < len; i++) {
		testArr3[i] = testArr1[i] * testArr2[i];
	}
	
	STOP_TIMER
	printf("multiplication test took: %f sec\n", sec);
	
	//	printf("%li\n", (long)testArr3);
	printf("%f\n", testArr3[100]);
}

#pragma mark - main

int main(int argc, const char * argv[]) {
	testFunc5();

	return 0;
}
















