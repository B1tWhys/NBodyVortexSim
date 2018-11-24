#include "constants.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define checkVal if (value == NULL) {\
    fprintf(stderr, "config warning: value expected for %s\n", keyword);\
    continue;\
}

// default values for all of the constants prior to being set from the config file
int TEST_CASE = 0;
char DRAW_CONSOLE = 0;
char DRAW_PNG = 0;
char SAVE_RAWDATA = 0;
char DATA_OUT_FILEPATH[255] = "";
char INITFNAME[255] = "";
int INIT_TIME_STEP = 0;
int CONSOLE_W = 200;
int CONSOLE_H = 100;
int IMAGE_W = 1000;
int IMAGE_H = 1000;
int VORTEX_DRAW_SIZE_CONST = 4;
int TRACER_DRAW_SIZE_CONST = 1;
#if !(defined(DOMAIN_SIZE_X) || defined(DOMAIN_SIZE_Y))
int DOMAIN_SIZE_X = 64;
int DOMAIN_SIZE_Y = 64;
#endif
float TIMESTEP_CONST = .01;
int RENDER_NTH_STEP = 5;
#ifndef NUMBER_OF_STEPS
int NUMBER_OF_STEPS = 1000;
#endif
int NUM_TRACERS = 1;
int NUM_VORT_INIT = 64;
int FIRST_SEED = -1;
char VORTEX_LIFECYCLE = 1;
float VORTEX_INTENSITY_SIGMA = 0.21233045007200477 * 2 * M_PI;
float VORTEX_SPAWN_RATE = 2.56;
int VORTEX_MERGE_RADIUS = 1;
int THREADCOUNT = 8;

void importConstants(char *filename) {
    if (filename == NULL) {
        filename = "./config";
    }

    char buff[255];
    size_t lineLen;
    
    FILE *configFile = fopen(filename, "r");
    while (!feof(configFile)) {
        fgets(buff, 255, configFile);
        
        lineLen = strlen(buff);

        char *crPtr = memchr(buff, '\n', lineLen);
        if (crPtr) {
            *crPtr = 0;
            --lineLen;
        }

        char *comment;
        if ((comment = memchr(buff, '#', lineLen)) != NULL) {
            lineLen = comment - buff;
        }

        if (lineLen <= 1) continue; // skip lines only containing a carraige return and entirely blank lines
        
        char *keyword = strtok(buff, " ");
        char *value = strtok(NULL, " ");
        char *valueEnd = (value == NULL) ? NULL : &value[strlen(value)+1];

        // for each possible keyword, set the corresponding variable to true
        if (strcmp(keyword, "TEST_CASE") == 0) {
            TEST_CASE = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "DRAW_CONSOLE") == 0) {
            DRAW_CONSOLE = 1;
        } else if (strcmp(keyword, "DRAW_PNG") == 0) {
            DRAW_PNG = 1;
        } else if (strcmp(keyword, "SAVE_RAWDATA") == 0) {
            SAVE_RAWDATA = 1;
        } else if (strcmp(keyword, "DATA_OUT_FILEPATH") == 0) {
            memcpy(DATA_OUT_FILEPATH, value, strlen(value)+1);
        } else if (strcmp(keyword, "INITFNAME") == 0) {
            memcpy(INITFNAME, value, strlen(value)+1);
        } else if (strcmp(keyword, "INIT_TIME_STEP") == 0) {
            INIT_TIME_STEP = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "CONSOLE_H") == 0) {
            CONSOLE_H = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "CONSOLE_W") == 0) {
            CONSOLE_W = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "IMAGE_W") == 0) {
            IMAGE_W = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "IMAGE_H") == 0) {
            IMAGE_H = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "VORTEX_DRAW_SIZE_CONST") == 0) {
            VORTEX_DRAW_SIZE_CONST = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "TRACER_DRAW_SIZE_CONST") == 0) {
            TRACER_DRAW_SIZE_CONST = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "DOMAIN_SIZE_X") == 0) {
            DOMAIN_SIZE_X = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "DOMAIN_SIZE_Y") == 0) {
            DOMAIN_SIZE_Y = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "TIMESTEP_CONST") == 0) {
            TIMESTEP_CONST = strtof(value, NULL);
        } else if (strcmp(keyword, "RENDER_NTH_STEP") == 0) {
            RENDER_NTH_STEP = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "NUMBER_OF_STEPS") == 0) {
            NUMBER_OF_STEPS = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "NUM_TRACERS") == 0) {
            NUM_TRACERS = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "NUM_VORT_INIT") == 0) {
            NUM_VORT_INIT = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "FIRST_SEED") == 0) {
            FIRST_SEED = strtol(value, NULL, 10);
        } else if (strcmp(keyword, "VORTEX_LIFECYCLE") == 0) {
            VORTEX_LIFECYCLE = strtof(value, NULL);
        } else if (strcmp(keyword, "VORTEX_INTENSITY_SIGMA") == 0) {
            VORTEX_INTENSITY_SIGMA = strtof(value, NULL);
        } else if (strcmp(keyword, "VORTEX_SPAWN_RATE") == 0) {
            VORTEX_SPAWN_RATE = strtof(value, NULL);
        } else if (strcmp(keyword, "VORTEX_MERGE_RADIUS") == 0) {
            VORTEX_MERGE_RADIUS = strtof(value, NULL);
        } else if (strcmp(keyword, "THREADCOUNT") == 0) {
            THREADCOUNT = strtol(value, NULL, 10);
        } else {
            fprintf(stderr, "error: could not parse config file line:\n%s\n", buff);
        }
    }
}
