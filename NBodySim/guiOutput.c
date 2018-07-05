#include "guiOutput.h"
#include "constants.h"
#include "main.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cairo.h>
#include <string.h>


void genFName(char *strBuffer, int frameNum) {
	strcpy(strBuffer, "./outputImages/frame_");
	char buffer[9];
	sprintf(buffer, "%i", frameNum);
	
	for (int i = 9-(int)strlen(buffer); i >= 0; i--) {
		strcat(strBuffer, "0");
	}
	strcat(strBuffer, buffer);
	strcat(strBuffer, ".png");
}

void drawToConsole(struct Vortex *vorts, int numVorts, struct Tracer *tracers) {
	printf("\033[3J");
	int pxArray[CONSOLE_W][CONSOLE_H];

	// -2 means draw a tracer there, -1 means empty, >= 0 means a vort is there.

	for (int y = 0; y < CONSOLE_H; y++) {
		for (int x = 0; x < CONSOLE_W; x++) {
			pxArray[x][y] = -1;
		}
	}

	for (int i = 0; i < numVorts; i++) {
		int xPxCoord = (int)(CONSOLE_W-1)*vorts[i].position[0]/DOMAIN_SIZE_X;
		int yPxCoord = (int)(CONSOLE_H-1)*vorts[i].position[1]/DOMAIN_SIZE_Y;

		pxArray[xPxCoord][yPxCoord] = vorts[i].vIndex;
	}

	for (int i = 0; i < NUM_TRACERS; i++) {
		int xPxCoord = (int)(CONSOLE_W-1)*tracers[i].position[0]/DOMAIN_SIZE_X;
		int yPxCoord = (int)(CONSOLE_H-1)*tracers[i].position[1]/DOMAIN_SIZE_Y;

		pxArray[xPxCoord][yPxCoord] = -2;
	}

	int count = 0;
	for (int y = CONSOLE_H; y >= 0; y--) {
		for (int x = 0; x < CONSOLE_W; x++) {
			if (pxArray[x][y] >= 0) {
				count++;
				for (int j = 0; j < ceil(log10(pxArray[x][y])); j++) printf("\010");
				printf("\033[91m");
				printf("%i", pxArray[x][y]);
				printf("\033[0m");
			} else if (pxArray[x][y] == -2) {
				printf(".");
			// } else if (x == CONSOLE_W/2 && y == CONSOLE_H/2) {
			// 	printf("+");
			// } else if (x == CONSOLE_W/2) {
			// 	printf("|");
			// } else if (y == CONSOLE_H/2) {
			// 	printf("-");
			} else {
				printf(" ");
			}
		}
		if (y != CONSOLE_H) printf("\n");
	}
}

void drawToFile(struct Vortex *vorts, int numVorts, struct Tracer *tracers, char filename[]) {
	cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_RGB24, IMAGE_W, IMAGE_H);
	cairo_t *cr = cairo_create(surface);
	
	cairo_scale(cr, 1, -1);
	cairo_translate(cr, 0, -cairo_image_surface_get_height(surface));
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_paint(cr);
	
	for (int vortIndex = 0; vortIndex < numVorts; vortIndex++) {
		struct Vortex *vort = &vorts[vortIndex];
		
		double xPos = IMAGE_W * vort->position[0] / DOMAIN_SIZE_X;
		double yPos = IMAGE_H * vort->position[1] / DOMAIN_SIZE_Y;
		double rad = fabs(vort->intensity) * VORTEX_DRAW_SIZE_CONST;
		
		if (vort->intensity > 0) {
			cairo_set_source_rgb(cr, 1, 0, 0);
		} else {
			cairo_set_source_rgb(cr, 0, 0, 1);
		}
		
		cairo_new_sub_path(cr);
		cairo_arc(cr, xPos, yPos, rad, 0, 2 * M_PI);
		cairo_close_path(cr);
		cairo_fill(cr);
	}
	
	cairo_set_source_rgb(cr, 1, 0, 1);
	
	int i = 0;
	for (int tracerIndex = 0; tracerIndex < NUM_TRACERS; tracerIndex++) {
		struct Tracer *tracer = &tracers[tracerIndex];
		
		double xPos = IMAGE_W * tracer->position[0] / DOMAIN_SIZE_X;
		double yPos = IMAGE_H * tracer->position[1] / DOMAIN_SIZE_Y;
		double rad = 2;
		
		cairo_new_sub_path(cr);
		cairo_arc(cr, xPos, yPos, rad, 0, 2 * M_PI);
		cairo_close_path(cr);
		cairo_fill(cr);
		i++;
	}
	
	cairo_surface_write_to_png(surface, filename);
	
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}
