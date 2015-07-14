#include "sor.h"
#include "boundary_val.h"
#include <math.h>

void sor(
  double omg,
  double dx,
  double dy,
  double dp,
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
	for(i = 1; i <= imax; i++) {
		for(j = 1; j<=jmax; j++) {
			if(Flag[i][j] == C_F) {
				P[i][j] = (1.0-omg)*P[i][j]
						+ coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
			}
		}
	}

	/* Compute the residual */
	rloc = 0;
	for(i = 1; i <= imax; i++) {
		for(j = 1; j <= jmax; j++) {
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


	/* Set values for horizontal boundaries (Neumann BC) */
	for(i = 1; i <= imax; i++) {
		P[i][0] = P[i][1];
		P[i][jmax+1] = P[i][jmax];
	}
	
	/* Set values for vertical boundaries (Dirichlet BC if dp != 0, Neumann BC otherwise) 
		If dp != 0, we assign it to the left boundary and set the right boundary to 0 */
	if(dp != 0){
		for(j = 0; j <= jmax; j++) {
			P[0][j] = 2 * dp - P[1][j]; 
			P[imax+1][j] = -P[imax][j];
		} 
	} else{
		for(j = 1; j <= jmax; j++) {
			P[0][j] = P[1][j];
			P[imax+1][j] = P[imax][j];
		} 
	}
	
	/* Set pressures of obstacle cells */
	for(i = 1; i <= imax; i++) {
		for(j = 1; j<=jmax; j++) {
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

