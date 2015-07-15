#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

enum BC_type {
	BC_NO_SLIP = 1, 
	BC_FREE_SLIP = 2, 
	BC_OUTFLOW = 3
};
/* Flags for obstacle surfaces */
enum B_xy {
	B_N = 64+1, 
	B_S = 64+2, 
	B_W = 64+4, 
	B_NW = 64+5, 
	B_SW = 64+6, 
	B_O = 64+8, 
	B_NO = 64+9, 
	B_SO = 64+10
};
/* Flags for free surface cells */
enum FS_xy {
	FS_N 	= 32+16+1,		/* Empty cell at north */
	FS_S	= 32+16+2,		/* Empty cell at south */
	FS_O	= 32+16+4,		/* Empty cell at east */
	FS_W	= 32+16+8		/* Empty cell at west */
};
/* Flags to indicate cell type */
enum C_type {
	C_F 	= 16,		/* Interior fluid cell */
	C_FS	= 48,		/* Free surface cell */
	C_E		= 0			/* Empty cell */
};

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
  int **flagField,
  double dx,
  double dy
);

void noslip(int imax, int jmax, double **U, double **V, int side);
void freeslip(int imax, int jmax, double **U, double **V, int side);
void outflow(int imax, int jmax, double **U, double **V, int side);
void movingwall(int imax, int jmax, double **U, double **V, int side);
void spec_boundary_val (char *problem, int imax, int jmax, double **U, double **V, double **P, double Re, double xlength, double ylength, double dp);

#endif
