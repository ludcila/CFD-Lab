#include "uvp.h"
#include <math.h>

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
	/* test_1 modified by SanYu */
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
	for(i = 1; i < imax; i++) {
		for(j = 1; j < jmax; j++) {			/* TODO: check i and j ranges */
			F[i][j] = U[i][j] + dt * ( 
				+ 1/Re*((U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx*dx)) 
				- 1/dx * ( pow((U[i][j]+U[i+1][j])/2, 2) - pow((U[i-1][j]+U[i][j])/2, 2) ) 
				- 1/dy * ( (V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])/4 - (V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])/4 )
				+ alpha/dy * ( abs(V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])/4 - abs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])/4 )
			);
			G[i][j] = V[i][j] + dt * ( 
				+ 1/Re*((V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy)) 
				- 1/dy * ( pow((V[i][j]+V[i][j+1])/2, 2) - pow((V[i][j-1]+V[i][j])/2, 2) ) 
				- 1/dx * ( (U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])/4 - (U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])/4 )
				- alpha/dx * ( abs(U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])/4 - abs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])/4 )
			);
		}
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
	for(i = 0; i < imax; i++) {
		for(j = 0; j < jmax; j++) {			/* TODO: check i and j ranges */
			RS[i][j] = 1/dt * ((F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy);
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
}
