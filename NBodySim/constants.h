//
//  constants.h
//  NBodySim
//


// NOTE: these constants can all be set in the config file. There are default
// values for them in constants.c which are used for any constants which aren't
// set in the config file

#ifndef constants_h
#define constants_h

#include <math.h>

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

extern int TEST_CASE;

// render mode selection
extern char DRAW_CONSOLE;
extern char DRAW_PNG;
extern char SAVE_RAWDATA;
extern char SAVE_RK_STEPS; // save each runge-kutta timestep, in addition to every normal timestep

// file names for the file to initialize the simulation from, and the filepath to
// write to. To disable initilzing from a source file, set INITFNAME to "".
extern char DATA_OUT_FILEPATH[255];
extern char INITFNAME[255];
extern int INIT_TIME_STEP;

extern int CONSOLE_W; // character dimensions to draw to console
extern int CONSOLE_H;

extern int IMAGE_W; // size of output frames
extern int IMAGE_H;
extern int VORTEX_DRAW_SIZE_CONST;
extern int TRACER_DRAW_SIZE_CONST;

extern int DOMAIN_SIZE_X; // size of one box of the simulation in units
extern int DOMAIN_SIZE_Y;

extern float TIMESTEP_CONST;
extern int RENDER_NTH_STEP; // speeds up the simulation display
#ifndef NUMBER_OF_STEPS
extern int NUMBER_OF_STEPS; // number of time steps to simulate. 0 to loop forever
#endif

extern int NUM_TRACERS; // NOTE: must be a square number
extern int NUM_VORT_INIT;
extern int FIRST_SEED; // seed the sim. -1 to use current unix time stamp

extern char VORTEX_LIFECYCLE; // controls whether vortices are merged/spawned
extern float VORTEX_INTENSITY_SIGMA;
extern float VORTEX_SPAWN_RATE;
extern int VORTEX_MERGE_RADIUS;

extern int nextVortID;

extern int THREADCOUNT;

void importConstants(char *);
#endif
