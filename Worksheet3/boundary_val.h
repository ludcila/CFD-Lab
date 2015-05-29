#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

/* Boundary condition types */

#define		BC_NO_SLIP		1
#define		BC_FREE_SLIP	2
#define		BC_OUTFLOW		3

/**
 * The boundary values of the problem are set.
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
);

void noslip(int imax, int jmax, double **U, double **V, int side);
void freeslip(int imax, int jmax, double **U, double **V, int side);
void outflow(int imax, int jmax, double **U, double **V, int side);
void movingwall(int imax, int jmax, double **U, double **V, int side);
void spec_boundary_val (char *problem, int imax, int jmax, double **U, double **V);

#endif
