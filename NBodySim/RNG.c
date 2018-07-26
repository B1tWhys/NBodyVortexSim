//
//  RNG.c
//  NBodySim
//
//  Created by Skyler Arnold on 7/24/18.
//

#include "RNG.h"
#include <stdlib.h>
#include <math.h>

/**
 generates a uniformly random double value in a range
 
 @param lowerBound the lower bound for the random number
 @param upperBound the upper bound for the random number
 
 @return returns a random uniform double in the given range
 */
double generateUniformRandInRange(double lowerBound, double upperBound) { // range is inclusive on both ends
	// TODO: investigate whether there's an off by 1 error in here
	double range = upperBound-lowerBound;
	double result = ((double)rand()/RAND_MAX)*range + lowerBound;
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





