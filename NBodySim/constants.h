//
//  constants.h
//  NBodySim
//
//  Created by Skyler Arnold on 6/15/18.
//  Copyright © 2018 Skyler Arnold. All rights reserved.
//

#ifndef constants_h
#define constants_h

/*
 -DNUMBER_OF_STEPS=0
 -DXCODE

 num:	test case				(val for NUM_VORT_INIT)
 0:		normal mode				any integer
 1:		two co-orbiting vorts	(2)
 2:		two parallel vorts		(2)
 3:		triangle of vorts		(3)
 */
#define TEST_CASE 1

#if !defined(CONSOLE_H) || !defined(CONSOLE_W)
#define CONSOLE_W 100 // character dimensions to draw to console
#define CONSOLE_H 50
#endif

#define WINDOW_W 1000 // px dimentions for output window (not implemented yet)
#define WINDOW_H 1000
#define IMAGE_W 1000 // size of output images (also not implemented yet)
#define IMAGE_H 1000

#define DOMAIN_SIZE_X 64 // size of one box of the simulation in units
#define DOMAIN_SIZE_Y 64

#define TIMESTEP .1
#define RENDER_NTH_STEP 30 // speeds up the simulation display
#ifndef NUMBER_OF_STEPS // don't change this line
#define NUMBER_OF_STEPS 0 // number of time steps to simulate. 10 to loop forever
#endif // ignore this line too

#define NUM_TRACERS 0 // works more neatly if this is a square number
#define NUM_VORT_INIT 2
#define FIRST_SEED -1 // seed the sim. -1 to use current unix time stamp
#define VORTEX_INTENSITY_INIT_UPPER_BOUND 20

#endif /* constants_h */
