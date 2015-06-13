#include "sor.h"
#include "boundary_val.h"
#include <math.h>

void sor(
  double omg,
  double dx,
  double dy,
  double dp,
  int il, int ir,
  int jb, int jt,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  int **Flag
) {

  int i,j;
  int count_num_fluid_cell = 0;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

	/* SOR iteration */
	for(i = il; i <= ir; i++) {
		for(j = jb; j <= jt; j++) {
			if(Flag[i][j] == C_F) {
				P[i][j] = (1.0-omg)*P[i][j]
						+ coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
			}
		}
	}

	/* Compute the residual */
	rloc = 0;
	for(i = il; i <= ir; i++) {
		for(j = jb; j <= jt; j++) {
			if(Flag[i][j] == C_F) {
				rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
						( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
				count_num_fluid_cell++;
			}
		}
	}
	rloc = rloc / count_num_fluid_cell;
	rloc = sqrt(rloc);
	
	/* Set residual */
	*res = rloc;

	
	/* BC for lower wall */
	if(jb == 0) {
		for(i = il; i <= ir; i++) {
			P[i][0] = P[i][1];
		}
	}

	/* BC for upper wall */
	if(jt == jmax) {
		for(i = il; i <= ir; i++) {
			P[i][jmax+1] = P[i][jmax];
		}
	}
	
	/* Set values for vertical boundaries (Dirichlet BC if dp != 0, Neumann BC otherwise) 
		If dp != 0, we assign it to the left boundary and set the right boundary to 0 */
	if(dp != 0){
		if(il == 0) {
			for(j = jb; j <= jt; j++) {
				P[0][j] = 2 * dp - P[1][j]; 
			} 
		}
		if(ir == imax) {
			for(j = jb; j <= jt; j++) {
				P[imax+1][j] = -P[imax][j];
			} 
		}
	} else{
		if(il == 0) {
			for(j = jb; j <= jt; j++) {
				P[0][j] = P[1][j];
			}
		}
		if(ir == imax) {
			for(j = jb; j <= jt; j++) {
				P[imax+1][j] = P[imax][j];
			}
		}
	}
	
	/* Set pressures of obstacle cells */
	for(i = il; i <= ir; i++) {
		for(j = jb; j <= jt; j++) {
			switch(Flag[i][j]) {
				case B_N:
					P[i][j] = P[i][j+1];
					break;
				case B_S:
					P[i][j] = P[i][j-1];
					break;
				case B_O:
					P[i][j] = P[i+1][j];
					break;
				case B_W:
					P[i][j] = P[i-1][j];
					break;
				case B_NW:
					P[i][j] = 0.5 * (P[i-1][j] + P[i][j+1]);
					break;
				case B_NO:
					P[i][j] = 0.5 * (P[i+1][j] + P[i][j+1]);
					break;
				case B_SW:
					P[i][j] = 0.5 * (P[i-1][j] + P[i][j-1]);
					break;
				case B_SO:
					P[i][j] = 0.5 * (P[i+1][j] + P[i][j-1]);
					break;
			}
		}
	}

}

