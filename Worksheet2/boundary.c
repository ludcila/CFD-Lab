#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdio.h>

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
	/* TODO Jae */
	int x, y, z, i;
	int cellIdx, fluidCellIdx;
	double *boundaryCell, *fluidCell;
	double densityOfFluidCell;
	double c_dot_u;
	int xx, yy, zz;
	int numGridPoints = xlength + 2;

	for(z = 0; z < numGridPoints; z++) {
		for(y = 0; y < numGridPoints; y++) {
			for(x = 0; x < numGridPoints; x++) {
			
				cellIdx = z * numGridPoints * numGridPoints + y * numGridPoints + x;
				
				/* Set boundary conditions if it is not a fluid cell */
				if(flagField[cellIdx] != 0) {
				
					boundaryCell = collideField + Q * cellIdx;
					
					for(i = 0; i < Q; i++) {
					
						/* Determine neighbor for this velocity direction */
						xx = x + LATTICEVELOCITIES[i][0];
						yy = y + LATTICEVELOCITIES[i][1];
						zz = z + LATTICEVELOCITIES[i][2];
						
						/* Check if neighbor is inside the domain and is not in the boundary */
						if(xx >= 1 && yy >= 1 && zz >= 1 && xx <= xlength && yy <= xlength && zz <= xlength) {
						
							fluidCellIdx = Q * (zz * numGridPoints * numGridPoints + yy * numGridPoints + xx);
							fluidCell = collideField + fluidCellIdx;
							
							if(flagField[cellIdx] == 1) {
								boundaryCell[i] = fluidCell[Q - i - 1];
							} else if(flagField[cellIdx] == 2) {
								computeDensity(fluidCell, &densityOfFluidCell);
								c_dot_u = LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2];
								boundaryCell[i] = fluidCell[Q - i - 1] + 2 * LATTICEWEIGHTS[i] * densityOfFluidCell * c_dot_u / (C_S * C_S);
							}
							
							/*if(x==1 && y==1 && z==xlength+1) {
								printf("(%d, %d, %d) ", LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]);
								printf("(%d, %d, %d) ", xx, yy, zz);
								printf("%f ", boundaryCell[i]);
							}*/
							
						} 
						
						
					}
				}
				
			}
		}
	}

}

