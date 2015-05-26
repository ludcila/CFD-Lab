#include "streaming.h"
#include "LBDefinitions.h"
#include <stdlib.h>
#include <stdio.h>

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
	
	/* Index for i-th direction */
	int i;

	/* Coordinates of cell to be streamed to */
  	int x, y, z;
  	
  	/* Coordinates of neighbor cells */
  	int xx, yy, zz;
  	
  	/* Pointer and array indices for current and neighbor cell */
	double *currentCell, *neighborCell;
	int cellIdx, neighborCellIdx;
	
	/* Total number of grid points per side */
	int numGridPoints = xlength + 2;
	
	/* Loop through all fluid cells */
	for(z=1;z<=xlength;z++){
		for(y=1;y<=xlength;y++){
			for(x=1;x<=xlength;x++){
			
				/* Determine current cell at (x, y, z) */
				cellIdx = Q * (z * numGridPoints * numGridPoints + y * numGridPoints + x);
				currentCell = streamField + cellIdx;
				
				for(i=0;i<Q;i++){
				
					/* Determine neighbor in direction i, at (xx, yy, zz) */
					xx = x - LATTICEVELOCITIES[i][0];
					yy = y - LATTICEVELOCITIES[i][1];
					zz = z - LATTICEVELOCITIES[i][2];
					neighborCellIdx = Q * (zz * numGridPoints * numGridPoints + yy * numGridPoints + xx);
					neighborCell = collideField + neighborCellIdx;
					
					/* Stream from this neighbor*/
					currentCell[i] = neighborCell[i];
					
				}
				
			}	
		}
	}
}

