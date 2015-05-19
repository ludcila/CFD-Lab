#ifndef _MAIN_C_
#define _MAIN_C_

#include <math.h>
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"
#include "helper.h"

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
	
	/* Read parameters from the input file */
	readParameters( szFileName,
		&xlength, &tau, &velocityWallX, &timesteps, &timestepsPerPlotting, argc, argv
	);
	velocityWall[0] = velocityWallX;
	
	/* Allocate memory */
	numCells = pow(xlength + 2, D);
	collideField = calloc(1, Q * numCells * sizeof(double));
	streamField = calloc(1, Q * numCells * sizeof(double)); 
	flagField = calloc(1, numCells * sizeof(int));
	
	/* Set initial values of the fields */
	initialiseFields(collideField, streamField, flagField, xlength);
	
	for(t = 0; t < timesteps; t++) {
	
		struct timespec start_t, end_t;
	
		/* Streaming step */
		double *swap = NULL;
		
		clock_gettime(CLOCK_MONOTONIC, &start_t);
		doStreaming(collideField, streamField, flagField, xlength);
		swap = collideField;
		collideField = streamField;
		streamField = swap;
		clock_gettime(CLOCK_MONOTONIC, &end_t);
		printf("t=%d\tStreaming\t(%f s)\n", t, elapsedTime(start_t, end_t));
		
		/* Collide step */
		clock_gettime(CLOCK_MONOTONIC, &start_t);
		doCollision(collideField, flagField, &tau, xlength);
		clock_gettime(CLOCK_MONOTONIC, &end_t);
		printf("t=%d\tCollision\t(%f s)\n", t, elapsedTime(start_t, end_t));
		
		/* Set boundary values */
		clock_gettime(CLOCK_MONOTONIC, &start_t);
		treatBoundary(collideField, flagField, velocityWall, xlength);
		clock_gettime(CLOCK_MONOTONIC, &end_t);
		printf("t=%d\tBoundary\t(%f s)\n", t, elapsedTime(start_t, end_t));

		/* Output */
		if(t % timestepsPerPlotting == 0) {
			clock_gettime(CLOCK_MONOTONIC, &start_t);
			writeVtkOutput(collideField, flagField, "output", t, xlength);
			clock_gettime(CLOCK_MONOTONIC, &end_t);
			printf("t=%d\tWrite output\t(%f s)\n", t, elapsedTime(start_t, end_t));
		}
		
		printf("----------------------------------------\n");
		
	}
	
	/* Free allocated memory */
	free(collideField);
	free(streamField);
	free(flagField);

	return 0;
}

#endif

