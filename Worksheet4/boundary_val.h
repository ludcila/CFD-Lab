#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

enum BC_type {
	BC_NO_SLIP = 1, 
	BC_FREE_SLIP = 2, 
	BC_OUTFLOW = 3
};
enum B_xy {
	B_N = 1, 
	B_S = 2, 
	B_W = 4, 
	B_NW = 5, 
	B_SW = 6, 
	B_O = 8, 
	B_NO = 9, 
	B_SO = 10
};
enum C_type {
	C_F = 16, 
	C_B = 0
};

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int il, int ir,
  int jb, int jt,
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

void spec_boundary_val (char *problem, int il, int ir, int jb, int jt, int imax, int jmax, double **U, double **V, double **P, double Re, double xlength, double ylength, double dp);

#endif
