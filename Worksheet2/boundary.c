#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdio.h>

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
	
	/* Index for i-th direction */
	int i;
	
	/* Boundary cell, coordinates and array index */
	double *boundaryCell;
	int x, y, z;
	int cellIdx;
	
	/* Fluid cell (neighbor of boundary cell), coordinates and array index */
	double *fluidCell;
	int xx, yy, zz;
	int fluidCellIdx;
	double densityOfFluidCell;
	double c_dot_u;
	
	/* Total number of grid points per side */
	int numGridPoints = xlength + 2;
	
	/* Loop through all cells (Does not look very efficient, but it is currently not our bottleneck) */
	for(z = 0; z < numGridPoints; z++) {
		for(y = 0; y < numGridPoints; y++) {
			for(x = 0; x < numGridPoints; x++) {
				
				/* Compute array index based on (x, y, z) coordinates */
				cellIdx = z * numGridPoints * numGridPoints + y * numGridPoints + x;
				
				/* Set boundary conditions if it is a boundary cell (i.e., not a fluid cell) */
				if(flagField[cellIdx] != CELL_TYPE_FLUID) {
				
					/* Get the boundary cell */
					boundaryCell = collideField + Q * cellIdx;
					
					/* Consider possible neighbors from all directions */
					for(i = 0; i < Q; i++) {
					
						/* Determine coordinates of neighbor at the i-th direction */
						xx = x + LATTICEVELOCITIES[i][0];
						yy = y + LATTICEVELOCITIES[i][1];
						zz = z + LATTICEVELOCITIES[i][2];
						
						/* Check whether the coordinates correspond to a fluid cell (x, y and z within [1, xlength]) */
						if(xx >= 1 && yy >= 1 && zz >= 1 && xx <= xlength && yy <= xlength && zz <= xlength) {
							
							/* Get the fluid cell */
							fluidCellIdx = Q * (zz * numGridPoints * numGridPoints + yy * numGridPoints + xx);
							fluidCell = collideField + fluidCellIdx;
							
							/* Treat no-slip boundary condition (Q-i-1 gives the index of the opposite velocity direction)*/
							if(flagField[cellIdx] == CELL_TYPE_BOUNDARY_NO_SLIP) {
								boundaryCell[i] = fluidCell[Q - i - 1];
								
							/* Treat moving wall boundary condition */
							} else if(flagField[cellIdx] == CELL_TYPE_BOUNDARY_MOVING_WALL) {
								computeDensity(fluidCell, &densityOfFluidCell);
								c_dot_u = LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2];
								boundaryCell[i] = fluidCell[Q - i - 1] + 2 * LATTICEWEIGHTS[i] * densityOfFluidCell * c_dot_u / (C_S * C_S);
							}
							
						} 
						
					}
				}
				
			}
		}
	}

}

