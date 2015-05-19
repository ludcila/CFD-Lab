#ifndef _MAIN_C_
#define _MAIN_C_

#include <math.h>
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"

int main(int argc, char *argv[]){
	const char *szFileName = "cavity100.dat";	

	double *collideField = NULL;
	double *streamField = NULL;
	int *flagField = NULL;
	
	int xlength;
	double tau;
	double velocityWallX, velocityWall[3] = {0};
	int timesteps;
	int timestepsPerPlotting;
	int t;
	int numCells;
	
	readParameters( szFileName,
		&xlength, &tau, &velocityWallX, &timesteps, &timestepsPerPlotting, argc, argv
	);/*add "const char *szFileName" by sanyu, timesteps setting 1000 at first, and timestepsPerPlotting 100*/
	velocityWall[0] = velocityWallX;
	
	/* Allocate memory */
	numCells = pow(xlength + 2, D);
	collideField = calloc(1, Q * numCells * sizeof(double));
	streamField = calloc(1, Q * numCells * sizeof(double)); 
	flagField = calloc(1, numCells * sizeof(int));
	
	initialiseFields(collideField, streamField, flagField, xlength);
	
	for(t = 0; t < timesteps; t++) {
	
		struct timespec start_t, end_t;
	
		/* Streaming step */
		double *swap = NULL;
		
		/*int x, y, z, i;*/
		clock_gettime(CLOCK_MONOTONIC, &start_t);
		doStreaming(collideField, streamField, flagField, xlength);
		swap = collideField;
		collideField = streamField;
		streamField = swap;
		clock_gettime(CLOCK_MONOTONIC, &end_t);
		printf("%d Streaming: %f\n", t, (double)end_t.tv_sec + 1e-9 * end_t.tv_nsec - start_t.tv_sec - 1e-9 * start_t.tv_nsec);
		
		/* Collide step */
		clock_gettime(CLOCK_MONOTONIC, &start_t);
		doCollision(collideField, flagField, &tau, xlength);
		clock_gettime(CLOCK_MONOTONIC, &end_t);
		printf("%d Collision: %f\n", t, (double)end_t.tv_sec + 1e-9 * end_t.tv_nsec - start_t.tv_sec - 1e-9 * start_t.tv_nsec);
		
		/* Set boundary values */
		clock_gettime(CLOCK_MONOTONIC, &start_t);
		treatBoundary(collideField, flagField, velocityWall, xlength);
		clock_gettime(CLOCK_MONOTONIC, &end_t);
		printf("%d Treat boundary: %f\n", t, (double)end_t.tv_sec + 1e-9 * end_t.tv_nsec - start_t.tv_sec - 1e-9 * start_t.tv_nsec);
		/*
		for(z=1;z<=xlength;z++){
			for(y=1;y<=xlength;y++){
				for(x=1;x<=xlength;x++){
					double density; int numGridPoints = xlength + 2;
					printf("(%d %d %d) ", x, y, z);
					computeDensity(collideField + Q*(z*numGridPoints*numGridPoints + y*numGridPoints + x), &density);
					printf("%f ", density);
					for(i = 0; i < Q; i++) {
						printf("%.3f ", collideField[Q*(z*numGridPoints*numGridPoints + y*numGridPoints + x) + i]);
					}
					if(fabs(density - 1) > 1e-3) {
						printf(" *");
					}
					printf("\n");
				}
			}
		}*/
		/* Output */
		if(t % timestepsPerPlotting == 0) {
			writeVtkOutput(collideField, flagField, "output", t, xlength);
		}
		
	}
	
	/* Free allocated memory */
	free(collideField);
	free(streamField);
	free(flagField);

	return 0;
}

#endif

