#include "boundary_val.h"

void boundaryvalues(
	int imax,
	int jmax,
	double **U,
	double **V
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

}
