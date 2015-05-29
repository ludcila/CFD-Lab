#include "sor.h"
#include <math.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  int **Flag
) {
  int i,j;
  int count_num_fluid_cell;		/* count neighboring fluid cells */
  double P_N, P_S, P_W, P_O;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax; j++) {
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }

  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <= jmax; j++) {
      rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }
  rloc = rloc/(imax*jmax);
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;


  /* set boundary values */
  for(i = 1; i <= imax; i++) {
    P[i][0] = P[i][1];
    P[i][jmax+1] = P[i][jmax];
  }
  for(j = 1; j <= jmax; j++) {
    P[0][j] = P[1][j];
    P[imax+1][j] = P[imax][j];
  }
	
  /* compute pressures of obstacle cells */
	for(i = 1; i <= imax; i++) {
		for(j = 1; j<=jmax; j++) {
			if(Flag[i][j] && 1){ /* B_N, north cell is fluid */
				P_N = P[i][j+1];
				count_num_fluid_cell++;
			} 

			if(Flag[i][j] && 2){ /* B_S, south cell is fluid */
				P_N = P[i][j-1];
				count_num_fluid_cell++;
			}

			if(Flag[i][j] && 4){ /* B_W, west cell is fluid */
				P_N = P[i-1][j];
				count_num_fluid_cell++;
			}

			if(Flag[i][j] && 8){ /* B_O, east cell is fluid */
				P_N = P[i+1][j];
				count_num_fluid_cell++;
			}
		}

		P[i][j] = (P_N + P_S + P_W + P_O) / count_num_fluid_cell;

		P_N = 0;
		P_S = 0;
		P_W = 0;
		P_O = 0;
		count_num_fluid_cell = 0;
	}

}

