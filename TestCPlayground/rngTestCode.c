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
//#include "RNG.h"
//#include "constants.h"

#include "ulcg.h"
#include "unif01.h"
#include "bbattery.h"

// old uniform RNG

double generateUniformRandInRange(double lowerBound, double upperBound) { // range is inclusive on both ends
	// TODO: investigate whether there's an off by 1 error in here
	double range = upperBound-lowerBound;
	double result = ((double)rand()/RAND_MAX)*range + lowerBound;
	return result;
}

/*
 this is a linear congruential pseudo random number generator. Obligatory wiki link: https://en.wikipedia.org/wiki/Linear_congruential_generator
 */
const long mod = 2147483647;
const long c = 0;
const long a = 16807;
//long lastX = (FIRST_SEED == -1) ? rand() : FIRST_SEED;
long lastX = 12345;
int randVal() {
	lastX = (a * lastX + c) % mod;
	return (unsigned int)(lastX & mod);
//	return (int)(lastX);
}

double randDouble() {
	return (double)randVal()/mod;
}

double oldRandDouble() {
	return generateUniformRandInRange(0, 1);
}

int main(int argc, const char * argv[]) {
	if (1) {
		unif01_Gen *gen = unif01_CreateExternGen01("newgen", randDouble);
		bbattery_SmallCrush(gen);
//		unif01_Gen *gen = unif01_CreateExternGen01("oldgen", oldRandDouble);
//		bbattery_SmallCrush(gen);
		ulcg_DeleteGen(gen);
	} else { for (int i = 0; i < 10; i ++) {
			printf("%lf\n", randDouble());
		}
	}
	return 0;
}
















