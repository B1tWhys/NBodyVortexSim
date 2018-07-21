//
//  TestCaseInitializers.h
//  NBodySim
//
//  Created by Skyler Arnold on 7/20/18.
//

#ifndef TestCaseInitializers_h
#define TestCaseInitializers_h

#include <stdio.h>
#include "constants.h"
#include "main.h"

void initialize_pair_orbit_test(struct Vortex *vortices, int n);
void initialize_pair_parallel_test(struct Vortex *vortices, int n);
void initialize_square_system(struct Vortex *vortices, int n);
void initialize_single_point(struct Vortex *vortices, int n);
void initialize_test_case_4(struct Vortex *vortices, int n);
void initialize_test(struct Vortex *vortices, int n);

#endif /* TestCaseInitializers_h */
