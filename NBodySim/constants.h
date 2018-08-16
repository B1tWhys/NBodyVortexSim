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
 2:		two parallel vorts	(2)
 3:		square of vorts		(4)
 4:		test case 4				(3)
 5:		single vort initializer	(1)
 6:		tracer accuracy test	(2)	(1 tracer)
 */
#define TEST_CASE 0

#if !defined(CONSOLE_H) || !defined(CONSOLE_W)
#define CONSOLE_W 200 // character dimensions to draw to console
#define CONSOLE_H 100
#endif

// render mode selection
// #define DRAW_CONSOLE
//#define DRAW_PDF
//#define SAVE_RAWDATA

//#define DATA_OUT_FILEPATH "./rawData"
#define INITFNAME "./rawData"
#define INIT_TIME_STEP 500

#define IMAGE_W 1000 // size of output frames
#define IMAGE_H 1000
#define VORTEX_DRAW_SIZE_CONST 4
#define TRACER_DRAW_SIZE_CONST 1

#if !(defined(DOMAIN_SIZE_X) || defined(DOMAIN_SIZE_Y))
#define DOMAIN_SIZE_X 64 // size of one box of the simulation in units
#define DOMAIN_SIZE_Y 64
#endif

#define TIMESTEP_CONST .01
#define RENDER_NTH_STEP 5 // speeds up the simulation display
#ifndef NUMBER_OF_STEPS
#define NUMBER_OF_STEPS 2 // number of time steps to simulate. 0 to loop forever
#endif

#define NUM_TRACERS 1 // NOTE: must be a square number
#define NUM_VORT_INIT 64
#define FIRST_SEED -1 // seed the sim. -1 to use current unix time stamp

#define VORTEX_LIFECYCLE
#define VORTEX_INTENSITY_SIGMA 0.21233045007200477 * 2 * M_PI
#define VORTEX_SPAWN_RATE 2.56
#define VORTEX_MERGE_RADIUS 1

#define THREADCOUNT 8

#endif
