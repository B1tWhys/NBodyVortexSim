//
//  RNG.c
//  NBodySim
//
//  Created by Skyler Arnold on 7/24/18.
//

#include "RNG.h"
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>

/*
 this is a linear congruential pseudo random number generator. Obligatory wiki link: https://en.wikipedia.org/wiki/Linear_congruential_generator
 values for mod, c and a taken from C++11's minstd_rand0
 */
const long mod = 2147483647; // 2^31 - 1
const long c = 0;
const long a = 16807;
long lastX;

long randVal() {
	lastX = (a * lastX + c) % mod;
	return (unsigned int)(lastX & mod);
}

/**
 generates a uniformly random double value in a range
 
 @param lowerBound the lower bound for the random number
 @param upperBound the upper bound for the random number
 
 @return returns a random uniform double in the given range
 */
double generateUniformRandInRange(double lowerBound, double upperBound) { // range is inclusive on both ends
	double range = upperBound-lowerBound;
	double result = ((double)randVal()/mod)*range + lowerBound;
	return result;
}

double nextNum = NAN;
double z1 = 0;
char generate = 0;
const double epsilon = DBL_EPSILON;

/**
 generates a random double from a normal distributuion centered at 0 with a given standard deviation
 
 @discussion this uses the Box-Muller transform to create pairs of random numbers.
 one is returned, and the other is stored for the next time that the function is called.
 
 @param sigma the standard deviation for the normal distribution
 */
double generateNormalRand(double sigma) {
	generate = !generate;
	
	if (!generate) {
		return z1 * sigma;
	}
	
	double u1, u2;
	do {
		u1 = generateUniformRandInRange(0, 1);
		u2 = generateUniformRandInRange(0, 1);
	} while (u1 <= DBL_EPSILON);
	
	double z0;
	z0 = sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(2 * M_PI * u2);
	return z0 * sigma;
}

/*
 algorithm Poisson generator based upon the inversion by sequential search:[37]
 init:
 	Let x ← 0, p ← e−λ, s ← p.
 	Generate uniform random number u in [0,1].
 while u > s do:
 	x ← x + 1.
 	p ← p * λ / x.
 	s ← s + p.
 return x.
 */

/**
 @discussion Computes a poisson distributed random integer numbers using the inverse transform sampling algorithm from wikipedia
 https://en.wikipedia.org/wiki/Poisson_distribution#Generating_Poisson-distributed_random_variables
 */
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





