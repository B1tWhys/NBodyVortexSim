#include "guiOutput.h"
#include "constants.h"
#include "main.h"
#include <stdio.h>
#include <math.h>
#include <cairo.h>


void drawToConsole(struct Vortex *vorts, int numVorts, struct Tracer *tracers, char filename[]) {
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
				printf("%i", pxArray[x][y]);
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

void drawToFile(struct Vortex *vorts, int numVorts, struct Tracer *tracers) {
	cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_RGB24, IMAGE_W, IMAGE_H);
	cairo_t *cr = cairo_create(surface);

	for (vortIndex = 0; vortIndex < numVorts; vortIndex++) {
		struct Vortex *vort = &vorts[vortIndex];

		double xPos = IMAGE_W*DOMAIN_SIZE_X / vort->position[0];
		double yPos = IMAGE_H*DOMAIN_SIZE_Y / vort->position[1];
		double rad = vort->intensity * VORTEX_DRAW_SIZE_CONST;

		cairo_new_sub_path(cr);
		cairo_arc(cr, xPos, yPos, radius, 0, M_2_PI);
		cairo_close_path(cr);
	}

	cairo_destroy(cr);
	cairo_surfae_write_to_png(surface, filename);
	cairo_surface_destroy(surface);
}
