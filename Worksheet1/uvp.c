#include "uvp.h"

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
	/* TODO by Luci */
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
	/* TODO by Luci */
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
