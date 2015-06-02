#include "boundary_val.h"
#include <string.h>
#include <stdio.h>

/* flags
* C_F: fluid cell
* C_B: obstacle cell
* B_xy: boundary cell
*/


/*
* wl, wr, wt, wb 
* 1 = noslip
* 2 = feeslip
* 3 = outflow
*/
void boundaryvalues(
	int imax,
	int jmax,
	double **U,
	double **V,
	int wl,
	int wr,
	int wt,
	int wb,
	int **Flag
) {

	int i, j;	
	
	/* Left wall */
	switch(wl) {
		case BC_NO_SLIP:
			for(j = 0; j <= jmax; j++) {
				U[0][j] = 0;
				V[0][j] = -V[1][j];
			}
			break;
		case BC_FREE_SLIP:
			for(j = 0; j <= jmax; j++) {
				U[0][j] = 0;
				V[0][j] = V[1][j];
			}
			break;
		case BC_OUTFLOW:
			for(j = 0; j <= jmax; j++) {
				U[0][j] = U[1][j];
				V[0][j] = V[1][j];
			}
			break;
	}

	/* Right wall */
	switch(wr) {
		case BC_NO_SLIP:
			for(j = 0; j <= jmax; j++) {
				U[imax][j] = 0;
				V[imax+1][j] = -V[imax][j];
			}
			break;
		case BC_FREE_SLIP:
			for(j = 0; j <= jmax; j++) {
				U[imax][j] = 0;
				V[imax+1][j] = V[imax][j];
			}
			break;
		case BC_OUTFLOW:
			for(j = 0; j <= jmax; j++) {
				U[imax][j] = U[imax-1][j];
				V[imax+1][j] = V[imax][j];
			}
			break;
	}
	
	/* Top wall */
	switch(wt) {
		case BC_NO_SLIP:
			for(i = 0; i <= imax; i++) {
				U[i][jmax+1] = -U[i][jmax];
				V[i][jmax] = 0;
			}
			break;
		case BC_FREE_SLIP:
			for(i = 0; i <= imax; i++) {
				U[i][jmax+1] = U[i][jmax];
				V[i][jmax] = 0;
			}
			break;
		case BC_OUTFLOW:
			for(i = 0; i <= imax; i++) {
				U[i][jmax+1] = U[i][jmax];
				V[i][jmax] = V[i][jmax-1];
			}
			break;
	}

	/* Bottom wall */
	switch(wb) {
		case BC_NO_SLIP:
			for(i = 0; i <= imax; i++) {
				U[i][0] = -U[i][1];
				V[i][0] = 0;
			}
			break;
		case BC_FREE_SLIP:
			for(i = 0; i <= imax; i++) {
				U[i][0] = U[i][1];
				V[i][0] = 0;
			}
			break;
		case BC_OUTFLOW:
			for(i = 0; i <= imax; i++) {
				U[i][0] = U[i][1];
				V[i][0] = V[i][1];
			}
			break;
	}
	
	/* Obstacles */
	for(i = 1; i < imax; i++) {
		for(j = 1; j < jmax; j++) {
		
			switch(Flag[i][j]) {
				case B_N:
					V[i][j] = 0;
					U[i][j] = -U[i][j+1];
					break;
				case B_S:
					V[i][j-1] = 0;
					U[i][j] = -U[i][j-1];
					break;
				case B_O:
					U[i][j] = 0;
					V[i][j] = -V[i+1][j];
					break;
				case B_W:
					U[i-1][j] = 0;
					V[i][j] = -V[i-1][j];
					break;
				case B_NW:
					V[i][j] = 0;
					U[i-1][j] = 0;
					U[i][j] = -U[i][j+1];
					break;
				case B_NO:
					U[i][j] = 0;
					V[i][j] = 0;
					break;
				case B_SW:
					U[i-1][j] = 0;
					V[i][j-1] = 0;
					U[i][j] = -U[i][j-1];
					V[i][j] = -V[i-1][j];
					break;
				case B_SO:
					U[i][j] = 0;
					V[i][j-1] = 0;
					V[i][j] = -V[i+1][j];
					break;
			}
			
		}
	}
	
}

void movingwall(int imax, int jmax, double **U, double **V, int side){

	int i, j;
	double U_wall = 1.0;
	if (side == 0){
		/* Left wall */
		for(j = 0; j <= jmax; j++) {
			U[0][j] = 0;
			V[0][j] = 2 * U_wall - V[1][j];
		}
	} else if (side == 1){
		/* Right wall */
		for(j = 0; j <= jmax; j++) {
			U[imax][j] = 0;
			V[imax+1][j] = 2 * U_wall - V[imax][j];
		}
	} else if (side == 2){
		/* Top */
		for(i = 0; i <= imax; i++) {
			U[i][jmax+1] = 2 * U_wall - U[i][jmax];
			V[i][jmax] = 0;
		}
	} else if (side == 3){
		/* Bottom */
		for(i = 0; i <= imax; i++) {
			U[i][0] = 2 * U_wall - U[i][1];
			V[i][0] = 0;
		}
	}
}

void spec_boundary_val (char *problem, int imax, int jmax, double **U, double **V, double **P, double Re, double xlength, double ylength, double dp){
	
	int j;
	if(dp == 0){
		if(strcmp(problem, "flow_over_step") == 0) {
			for(j = jmax/2 + 1; j <= jmax; j++) {
				U[0][j] = 1;
				V[0][j] = 0;
			}
		} else if(strcmp(problem, "karman_vortex_street") == 0) {
			for(j = 0; j <= jmax; j++) {
				U[0][j] = 1;
				V[0][j] = 0;
			}
		} else if(strcmp(problem, "plane_shear_flow") == 0) {
			for(j = 0; j <= jmax; j++) {
				U[0][j] = 1;
				V[0][j] = 0;
			}
		}
	}
}
