#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"


void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
	/* TODO Jae */
	int x, y, z, i;
	int cellIdx;
	double *currentCell;
	double density;
	double velocity[3];
	double c_dot_u;			

	for(z = 1; z <= xlength; z++) {
		for(y = 1; y <= xlength; y++) {
			for(x = 1; x <= xlength; x++) {
				cellIdx = Q * (z * xlength * xlength + y * xlength + x);  
				currentCell = collideField + cellIdx;
				
				computeDensity(currentCell, &density);

				if(*(flagField + cellIdx) == 1)
				for ( i = 0; i < Q; i++){
					c_dot_u = LATTICEVELOCITIES[i][0] * velocity[0] + LATTICEVELOCITIES[i][1] * velocity[1] + LATTICEVELOCITIES[i][2] * velocity[2];
					currentCell[i] = currentCell[i] + 2 * LATTICEWEIGHTS[i] * density * c_dot_u * wallVelocity[1] / (C_S * C_S);				
				}
				else if(*(flagField + cellIdx) == 2)
				for ( i = 0 ; i < Q; i++){
					currentCell[i] = currentCell[i] + 2 * LATTICEWEIGHTS[i] * density * c_dot_u * wallVelocity[2] / (C_S * C_S);
					computeDensity(currentCell, &density);
					currentCell[i] = currentCell[i] + 2 * LATTICEWEIGHTS[i] * density * c_dot_u * wallVelocity[2] / (C_S * C_S);
					
				}

			}
		}
	}	
}

