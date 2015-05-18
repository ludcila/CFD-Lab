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
	double velocityWall[3];
	int timesteps;
	int timestepsPerPlotting;
	int t;
	int numCells;
	
	readParameters( szFileName,
		&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv
	);/*add "const char *szFileName" by sanyu, timesteps setting 1000 at first, and timestepsPerPlotting 100*/
	
	/* Allocate memory */
	numCells = pow(xlength + 2, D);
	collideField = calloc(1, Q * numCells * sizeof(double));
	streamField = calloc(1, Q * numCells * sizeof(double)); 
	flagField = calloc(1, numCells * sizeof(int)); /*the dimension of flagField needs to be checked*/
	
	initialiseFields(collideField, streamField, flagField, xlength);
	
	for(t = 0; t < timesteps; t++) {
	
		/* Streaming step */
		double *swap = NULL;
		doStreaming(collideField, streamField, flagField, xlength);
		swap = collideField;
		collideField = streamField;
		streamField = swap;
		
		/* Collide step */
		doCollision(collideField, flagField, &tau, xlength);
		
		/* Set boundary values */
		treatBoundary(collideField, flagField, velocityWall, xlength);
		
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
