//
//  constants.h
//  NBodySim
//
//

#ifndef constants_h
#define constants_h

/*
 Test case numbering info:
 
 num:	test case				(val for NUM_VORT_INIT)
 0:		normal mode				(*)
 1:		two co-orbiting vorts	(2)
 2:		two parallel vorts		(2)
 3:		square of vorts			(3)
 4:		test case 4				(3)
 5:		single vort initializer	(1)
 6:		traer accuracy test		(2) (1 tracer)
 */
#define TEST_CASE 0

#if !defined(CONSOLE_H) || !defined(CONSOLE_W)
#define CONSOLE_W 212 // character dimensions to draw to console
#define CONSOLE_H 106
#endif

// render mode selection
//#define DRAW_CONSOLE
#define DRAW_PDF

#define WINDOW_W 1000 // px dimentions for output window (not implemented yet)
#define WINDOW_H 1000
#define IMAGE_W 500 // size of output images (also not implemented yet)
#define IMAGE_H 500
#define VORTEX_DRAW_SIZE_CONST 5
#define TRACER_DRAW_SIZE_CONST 2

#define DOMAIN_SIZE_X 64 // size of one box of the simulation in units
#define DOMAIN_SIZE_Y 64

#define TIMESTEP_CONST .5
#define RENDER_NTH_STEP 10 // speeds up the simulation display
#ifndef NUMBER_OF_STEPS // don't change this line
#define NUMBER_OF_STEPS 1000 // number of time steps to simulate. 0 to loop forever
#endif // ignore this line too

#define NUM_TRACERS 50176 // NOTE: must be a square number
#define NUM_VORT_INIT 150
#define FIRST_SEED -1 // seed the sim. -1 to use current unix time stamp
#define VORTEX_INTENSITY_SIGMA 0.21233045007200477
#define VORTEX_MERGE_RADIUS_CUTOFF .1

#define THREADCOUNT 8

#endif /*
constants_h

void deleteVortex(struct Vortex *vort, double *vortexRadii, int *numActiveVorts, struct Vortex *vorts, int numVorts) {
 */
