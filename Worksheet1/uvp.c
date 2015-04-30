#include "uvp.h"
#include <math.h>
#include <stdlib.h>
#include "helper.h"
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
	/* TODO by San Yu */

	int i, j;
	double u_max = U[1][1], v_max = V[1][1];
	double dt_min;

	for (i = 1; i<imax; i++){				/*double for-loop for max velocity in x&y direction*/
		for (j = 1; j<jmax; j++){
			if (U[i][j]>u_max)
				u_max = U[i][j];
			if (V[i][j]>v_max)
				v_max = V[i][j];

		}
	}

	dt_min = fmin(dy / v_max, dx / u_max);
	dt_min = fmin(Re*0.5 * 1 / (1 / (dx*dx) + 1 / (dy*dy)), dt_min);
	*dt = tau*dt_min;
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
	double **G
	) {
	int i, j;
	for(j = 1; j <= jmax; j++) {
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}
	for (i = 1; i < imax; i++) {
		for (j = 1; j <= jmax; j++) {			/* TODO: check i and j ranges */
			F[i][j] = U[i][j] + dt * (
				+ 1 / Re*((U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) / (dx*dx) + (U[i][j+1] - 2 * U[i][j] + U[i][j-1]) / (dy*dy))
				- 1 / dx * (pow((U[i][j] + U[i + 1][j]) / 2, 2) - pow((U[i - 1][j] + U[i][j]) / 2, 2))
				- alpha / dx * (abs(U[i][j] + U[i + 1][j])*(U[i][j] - U[i + 1][j]) / 4 - abs(U[i - 1][j] + U[i][j])*(U[i - 1][j] - U[i][j]) / 4)
				- 1 / dy * ((V[i][j] + V[i + 1][j])*(U[i][j] + U[i][j + 1]) / 4 - (V[i][j - 1] + V[i + 1][j - 1])*(U[i][j - 1] + U[i][j]) / 4)
				- alpha / dy * (abs(V[i][j] + V[i + 1][j])*(U[i][j] - U[i][j + 1]) / 4 - abs(V[i][j - 1] + V[i + 1][j - 1])*(U[i][j - 1] - U[i][j]) / 4)
			);
		}
	}
	for (i = 1; i <= imax; i++) {
		G[i][0] = V[i][0];
		for (j = 1; j < jmax; j++) {
			G[i][j] = V[i][j] + dt * (
				+ 1 / Re*((V[i][j + 1] - 2 * V[i][j] + V[i][j - 1]) / (dy*dy) + (V[i+1][j] - 2 * V[i][j] + V[i - 1][j]) / (dx*dx))
				- 1 / dy * (pow((V[i][j] + V[i][j + 1]) / 2, 2) - pow((V[i][j - 1] + V[i][j]) / 2, 2))
				- alpha / dy * (abs(V[i][j] + V[i][j + 1])*(V[i][j] - V[i][j+1]) / 4 - abs(V[i][j - 1] + V[i][j])*(V[i][j - 1] - V[i][j]) / 4)
				- 1 / dx * ((U[i][j] + U[i][j + 1])*(V[i][j] + V[i + 1][j]) / 4 - (U[i - 1][j] + U[i - 1][j + 1])*(V[i - 1][j] + V[i][j]) / 4)
				- alpha / dx * (abs(U[i][j] + U[i][j + 1])*(V[i][j] - V[i + 1][j]) / 4 - abs(U[i - 1][j] + U[i - 1][j + 1])*(V[i - 1][j] - V[i][j]) / 4)
			);
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
		for (j = 1; j <= jmax; j++) {			/* TODO: check i and j ranges */
			RS[i][j] = 1 / dt * ((F[i][j] - F[i - 1][j]) / dx + (G[i][j] - G[i][j - 1]) / dy);
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
	double **P
	) {



	/* TODO by Jae */
	unsigned int i, j;
	for (i = 1; i < imax-1; i++){				/*double for-loop for max velocity in x&y direction*/
		for (j = 1; j <= jmax; j++){
			U[i][j] = F[i][j] - dt * (P[i+1][j] - P[i][j]) / dx;
			
		}
	}
	for (i = 1; i <=imax; i++){				/*double for-loop for max velocity in x&y direction*/
		for (j = 1; j < jmax-1; j++){
			V[i][j] = G[i][j] - dt * (P[i][j+1] - P[i][j]) / dy;
		}
	}
}
