#include "uvp.h"
#include <math.h>
#include <stdlib.h>
#include "helper.h"

enum B_xy {B_NULL = 0, B_N = 1, B_S = 2, B_W = 4, B_NW = 5, B_SW = 6, B_NSW = 7, B_O = 8, B_NO = 9, B_SO = 10, B_NSO = 11, B_NWO = 13, B_SWO = 14};

void calculate_dt(
	double Re,
	double tau,
	double *dt,
	double dx,
	double dy,
	int imax,
	int jmax,
	double **U,
	double **V
) {

	int i, j;
	double u_max = U[1][1], v_max = V[1][1];
	double dt_min;

	for (i = 1; i <= imax; i++) {				/*double for-loop for max velocity in x&y direction*/
		for (j = 1; j <= jmax; j++) {
			if (fabs(U[i][j]) > u_max)
				u_max = fabs(U[i][j]);
			if (fabs(V[i][j]) > v_max)
				v_max = fabs(V[i][j]);
		}
	}

	dt_min = fmin(dy / v_max, dx / u_max);
	dt_min = fmin(Re * 0.5 * 1 / (1 / (dx * dx) + 1 / (dy * dy)), dt_min);
	*dt = tau * dt_min;
	
}


void calculate_fg(
	double Re,
	double GX,
	double GY,
	double alpha,
	double dt,
	double dx,
	double dy,
	int imax,
	int jmax,
	double **U,
	double **V,
	double **F,
	double **G,
	int **Flag
) {
	int i, j;
	for(j = 1; j <= jmax; j++) {
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}
	for (i = 1; i <= imax - 1; i++) {
		for (j = 1; j <= jmax; j++) {
			if(Flag[i][j] & B_W) { /*B_W, west cell is fluid*/
				F[i-1][j] = U[i-1][j];
			} else if(Flag[i][j] & B_O) { /*B_O, east cell is fluid*/
				F[i][j] = U[i][j];
			} else {
			F[i][j] = U[i][j] + dt * (
			              + 1 / Re * ((U[i+1][j] - 2 * U[i][j] + U[i-1][j]) / (dx * dx) + (U[i][j+1] - 2 * U[i][j] + U[i][j-1]) / (dy * dy))
			              - 1 / dx * (pow((U[i][j] + U[i+1][j]) / 2, 2) - pow((U[i-1][j] + U[i][j]) / 2, 2))
			              - alpha / dx * (fabs(U[i][j] + U[i+1][j]) * (U[i][j] - U[i+1][j]) / 4 - fabs(U[i-1][j] + U[i][j]) * (U[i-1][j] - U[i][j]) / 4)
			              - 1 / dy * ((V[i][j] + V[i+1][j]) * (U[i][j] + U[i][j+1]) / 4 - (V[i][j-1] + V[i+1][j-1]) * (U[i][j-1] + U[i][j]) / 4)
			              - alpha / dy * (fabs(V[i][j] + V[i+1][j]) * (U[i][j] - U[i][j+1]) / 4 - fabs(V[i][j-1] + V[i+1][j-1]) * (U[i][j-1] - U[i][j]) / 4)
			              + GX);
			}
		}
	}
	for (i = 1; i <= imax; i++) {
		G[i][0] = V[i][0];
		for (j = 1; j <= jmax - 1; j++) {
			if(Flag[i][j] & B_N) { /*B_N, north cell is fluid*/
				G[i][j] = V[i][j];
			} else if(Flag[i][j] & B_S) { /*B_S, south cell is fluid*/
				G[i][j-1] = V[i][j-1];
			} else {			

				G[i][j] = V[i][j] + dt * (
					      + 1 / Re * ((V[i][j+1] - 2 * V[i][j] + V[i][j-1]) / (dy * dy) + (V[i+1][j] - 2 * V[i][j] + V[i-1][j]) / (dx * dx))
					      - 1 / dy * (pow((V[i][j] + V[i][j+1]) / 2, 2) - pow((V[i][j-1] + V[i][j]) / 2, 2))
					      - alpha / dy * (fabs(V[i][j] + V[i][j+1]) * (V[i][j] - V[i][j+1]) / 4 - fabs(V[i][j-1] + V[i][j]) * (V[i][j-1] - V[i][j]) / 4)
					      - 1 / dx * ((U[i][j] + U[i][j+1]) * (V[i][j] + V[i+1][j]) / 4 - (U[i-1][j] + U[i-1][j+1]) * (V[i-1][j] + V[i][j]) / 4)
					      - alpha / dx * (fabs(U[i][j] + U[i][j+1]) * (V[i][j] - V[i+1][j]) / 4 - fabs(U[i-1][j] + U[i-1][j+1]) * (V[i-1][j] - V[i][j]) / 4)
					      + GY);
			}
		}
		G[i][jmax] = V[i][jmax];
	}
}


void calculate_rs(
	double dt,
	double dx,
	double dy,
	int imax,
	int jmax,
	double **F,
	double **G,
	double **RS
) {

	int i, j;
	for (i = 1; i <= imax; i++) {
		for (j = 1; j <= jmax; j++) {
			RS[i][j] = 1 / dt * ((F[i][j] - F[i-1][j]) / dx + (G[i][j] - G[i][j-1]) / dy);
		}
	}

}


void calculate_uv(
	double dt,
	double dx,
	double dy,
	int imax,
	int jmax,
	double **U,
	double **V,
	double **F,
	double **G,
	double **P,
	int **Flag
) {
/*
	// treat v
	for (i = 1; i <= imax - 1; i++) {
		for (j = 1; j <= jmax; j++) {
			if(Flag[i][j] && 1){//B_O, north cell is fluid
				V[i][j] = 0;
			} else if(!(Flag[i][j] && 1)){
				V[i][j] = -V[i-1][j];
			}

			if(Flag[i][j] && 2){//B_O, south cell is fluid
				V[i][j-1] = 0;	
			} else if(!(Flag[i][j] && 2)){
				V[i][j-1] = -V[i-1][j-1];
			}

			if(Flag[i][j] && 4){//B_O, west cell is fluid
				V[i][j] = -V[i+1][j];

			} else if(!(Flag[i][j] && 4)){
			}

			if(Flag[i][j] && 8){//B_O, east cell is fluid
				V[i][j] = -V[i+1][j];
			} else if(!(Flag[i][j] && 8)){
			}
			
		}
	}
	// treat u
	for (i = 1; i <= imax; i++) {
		for (j = 1; j <= jmax - 1; j++) {
			if(Flag[i][j] && 1){//B_O, north cell is fluid
			} else if(!(Flag[i][j] && 1)){
			}
				
			if(Flag[i][j] && 2){//B_O, south cell is fluid
			} else if(!(Flag[i][j] && 2)){
			}

			if(Flag[i][j] && 4){//B_O, west cell is fluid
				U[i-1][j] = 0;
				U[i-1][j] = -U[i-1][j+1];	
			} else if(!(Flag[i][j] && 4)){
			}

			if(Flag[i][j] && 8){//B_O, east cell is fluid
				U[i][j] = 0;
			} else if(!(Flag[i][j] && 8)){
				U[i][j] = -U[i][j+1];
			}
			
		}
	}
*/
	unsigned int i, j;
	for (i = 1; i <= imax - 1; i++) {
		for (j = 1; j <= jmax; j++) {
			U[i][j] = F[i][j] - dt * (P[i+1][j] - P[i][j]) / dx;

		}
	}
	for (i = 1; i <= imax; i++) {
		for (j = 1; j <= jmax - 1; j++) {
			V[i][j] = G[i][j] - dt * (P[i][j+1] - P[i][j]) / dy;
		}
	}
	
}
