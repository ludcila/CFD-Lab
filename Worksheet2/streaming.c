#include "streaming.h"
#include "LBDefinitions.h"
#include <stdlib.h>
#include <stdio.h>

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
  	int x,y,z;
  	int xx, yy, zz;
	int i;
	int cellIdx, neighborCellIdx;
	double *currentCell, *neighborCell;
	int numGridPoints = xlength + 2;

	for(z=1;z<=xlength;z++){
		for(y=1;y<=xlength;y++){
			for(x=1;x<=xlength;x++){
				double sum = 0;
				cellIdx = Q * (z * numGridPoints * numGridPoints + y * numGridPoints + x);
				currentCell = streamField + cellIdx;
				for(i=0;i<Q;i++){
					xx = x -LATTICEVELOCITIES[i][0];
					yy = y -LATTICEVELOCITIES[i][1];
					zz = z -LATTICEVELOCITIES[i][2];
					neighborCellIdx = Q * (zz * numGridPoints * numGridPoints + yy * numGridPoints + xx);
					neighborCell = collideField + neighborCellIdx;
					currentCell[i] = neighborCell[i];
					sum += currentCell[i];
				}
				/*printf("%f ", sum);*/
			}	
		}
	}
}

