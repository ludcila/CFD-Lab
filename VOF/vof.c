#include "helper.h"
#include "vof.h"

/* Set the initial fluid fraction based on the pgm file */
void init_fluidFraction(int **pgm, double **fluidFraction, int imax, int jmax) {
	int i, j;
	for(i = 1; i <= imax; i++) {
		for(j = 1; j <= jmax; j++) {
			if(pgm[i][j] == 1) {
				fluidFraction[i][j] = 1;
			}
		}
	}
}

/* Bookeeping adjustments as described by Hirt and Nichols */
void adjust_fluidFraction(double **fluidFraction, int **flagField, double epsilon, int imax, int jmax) {

	int i, j;
	int O_empty, W_empty, N_empty, S_empty;
	
	/* TODO: Take care of arbitrary geometry */
	
	for(i = 1; i <= imax; i++) {
		for(j = 1; j <= jmax; j++) {
		
			/* Treat empty cells */
			if(fluidFraction[i][j] < epsilon) {
				fluidFraction[i][j] = 0;
				flagField[i][j] = C_E;
			}
			
			/* Treat cells that contain fluid (either full or a fraction) */
			else {
			
				/* Initially set as full fluid cell (no empty neighbors) */
				flagField[i][j] = C_F;
				
				/* Determine empty neighbors; the min/max prevent it from checking neighbors on the domain boundary */
				O_empty = fluidFraction[min(i+1, imax)][j] < epsilon;
				W_empty = fluidFraction[max(i-1, 1)][j] < epsilon;
				N_empty = fluidFraction[i][min(j+1, jmax)] < epsilon;
				S_empty = fluidFraction[i][max(j-1, 1)] < epsilon;
				
				/* Toggle bit for each empty neighbor empty neighbors
					if there are no empty neighbors, the flag will remain unchanged (i.e., C_F)
					if there are empty neighbors, the flag will be a combination of FS_x flags
				*/
				flagField[i][j] = flagField[i][j] | (O_empty * FS_O) | (W_empty * FS_W) | (N_empty * FS_N) | (S_empty * FS_S);
				
				/* If it is a full fluid cell, we reduce the fluid fraction by 1.1*epsilon for each empty neighbor */
				if(fluidFraction[i][j] > 1-epsilon)
					fluidFraction[i][j] = 1 - 1.1 * epsilon * (O_empty + W_empty + N_empty + S_empty);
				
			}
			
		}
	}

}

/* Determine orientation of the free surface */
void calculate_freeSurfaceOrientation(double **fluidFraction, int **flagField, double **dFdx, double **dFdy, double dx, double dy, int imax, int jmax) {

	int i, j;
	
	/* Copy fluidFraction to boundary cells, just to make computations easier */
	for(i = 1; i <= imax; i++) {
		fluidFraction[i][0] = fluidFraction[i][1];
		fluidFraction[i][jmax+1] = fluidFraction[i][jmax];
	}
	for(j = 1; j <= jmax; j++) {
		fluidFraction[0][j] = fluidFraction[1][j];
		fluidFraction[imax+1][j] = fluidFraction[imax][j];
	}
	
	/* (Assuming uniform grid in each axis) */
	
	for(i = 1; i <= imax; i++) {
		for(j = 1; j <= jmax; j++) {
			if(flagField[i][j] & C_FS) {
				dFdx[i][j] = 0.5 * ((fluidFraction[i+1][j-1] + fluidFraction[i+1][j] + fluidFraction[i+1][j+1]) - (fluidFraction[i-1][j-1] + fluidFraction[i-1][j] + fluidFraction[i-1][j+1])) / dx;
				dFdy[i][j] = 0.5 * ((fluidFraction[i-1][j+1] + fluidFraction[i][j+1] + fluidFraction[i+1][j+1]) - (fluidFraction[i-1][j-1] + fluidFraction[i][j-1] + fluidFraction[i+1][j-1])) / dy;
			}
		}
	}
	
}

/* Timestepping for the fluid fraction field */
void calculate_fluidFraction(
	double **fluidFraction, 
	double **U, 
	double **V, 
	double **dFdx, 
	double **dFdy, 
	int imax,
	int jmax,
	double dx, 
	double dy, 
	double dt
) {

}
