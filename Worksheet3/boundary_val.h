#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

/* Boundary condition types */

#define		BC_NO_SLIP		1
#define		BC_FREE_SLIP	2
#define		BC_OUTFLOW		3

enum B_xy {B_NULL = 0, B_N = 1, B_S = 2, B_W = 4, B_NW = 5, B_SW = 6, B_NSW = 7, B_O = 8, B_NO = 9, B_SO = 10, B_NSO = 11, B_NWO = 13, B_SWO = 14};

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
void spec_boundary_val (char *problem, int imax, int jmax, double **U, double **V, double **P, double Re, double xlength, double ylength, double dp);

#endif
