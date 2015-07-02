#include "boundary_val.h"
#include "vof.h"
#include <stdio.h>

void boundaryvalues(
	int imax,
	int jmax,
	double **U,
	double **V,
	int **flagField,
	double dx,
	double dy
) {

	int i, j;
	double U_wall = 0.0;
	
	/* Left wall */
	for(j = 0; j <= jmax; j++) {
		U[0][j] = 0;
		V[0][j] = -V[1][j];
	}
	
	/* Right wall */
	for(j = 0; j <= jmax; j++) {
		U[imax][j] = 0;
		V[imax+1][j] = -V[imax][j];
	}
	
	/* Floor */
	for(i = 0; i <= imax; i++) {
		U[i][0] = -U[i][1];
		V[i][0] = 0;
	}
	
	/* Ceiling */
	for(i = 0; i <= imax; i++) {
		U[i][jmax+1] = 2 * U_wall - U[i][jmax];
		V[i][jmax] = 0;
	}
	
	/* Boundary conditions for free surfaces */
	
	for(i = 1; i <= imax; i++) {
		for(j = 1; j <= jmax; j++) {
			
			if((flagField[i][j] & C_FS) == C_FS) {
				
				
				/* Treat free surface cells with only one empty neighbor */
				
				/* Satisfy du/dx + du/dy = 0 */
				if(flagField[i][j] == FS_W) {
					U[i-1][j] = dx * U[i][j] + dx * (V[i][j] - V[i][j-1]) / dy;
				} else if(flagField[i][j] == FS_O) {
					U[i][j] = dx * U[i-1][j] - dx * (V[i][j] - V[i][j-1]) / dy;
				} else if(flagField[i][j] == FS_S) {
					V[i-1][j] = dy * V[i][j] + dy * (U[i][j] - U[i-1][j]) / dx;
				} else if(flagField[i][j] == FS_N) {
					V[i][j] = dy * V[i][j-1] - dy * (U[i][j] - U[i-1][j]) / dx;
				} else {
				
					/* Treat free surfaces with two or more empty neighbors */
			
					/* Satisfy du/dx = 0 */
					if(flagField[i][j] & FS_W) {
						U[i-1][j] = U[i][j];
					} else if(flagField[i][j] & FS_O) {
						U[i][j] = U[i-1][j];
					}
				
					/* Satisfy dv/dy = 0 */
					if(flagField[i][j] & FS_S) {
						V[i][j-1] = V[i][j];
					} else if(flagField[i][j] & FS_N) {
						V[i][j] = V[i][j-1];
					}
					
				}
			}
			
		}
	}
	

}
