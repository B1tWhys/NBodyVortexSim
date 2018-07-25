//
//  RNG.h
//  NBodySim
//
//  Created by Skyler Arnold on 7/24/18.
//

#ifndef RNG_h
#define RNG_h

#include <stdio.h>

double generateUniformRandInRange(double lowerBound, double upperBound);
double generateNormalRand(double sigma);
int generatePoissonRand(double k, double L, double x);

#endif /* RNG_h */
