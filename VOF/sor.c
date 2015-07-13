#include "sor.h"
#include "vof.h"
#include <math.h>
#include <stdio.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double **fluidFraction,
  int **flagField,
  double **dFdx,
  double **dFdy,
  double *res
) {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
  int numFluidCells = 0;
  double m, eta, d, dnorm, dpar;
  double Fa, Fb;
  double Pfluid;
  double Psurface = 0;

  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax; j++) {
    	if(flagField[i][j] == C_F) {
		  P[i][j] = (1.0-omg)*P[i][j]
		          + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
		  numFluidCells++;
        }
    }
  }

  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <= jmax; j++) {
    	if(flagField[i][j] == C_F) {
      		rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
        }
    }
  }
  rloc = rloc/numFluidCells;
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;
  
	for(i = 1; i <= imax; i++) {
		for(j = 1; j <= jmax; j++) {
		
			if((flagField[i][j] & C_FS) == C_FS) {
			
				/* Determine orientation and pressure of neighbor fluid cell based on orientation */
				if(fabs(dFdy[i][j]) > fabs(dFdx[i][j])) { 					/* More horizontal than vertical */
					dnorm = dy;
					dpar = dx;
					if((flagField[i][j] & FS_N) == FS_N) {		/* If north neighbor is empty, fluid pressure is from south neighbor */
						Pfluid = P[i][j-1];
					} else if((flagField[i][j] & FS_S) == FS_S) {	/* If south neighbor is empty, fluid pressure is from north neighbor */
						Pfluid = P[i][j+1];
					}
				} else { 							/* More vertical than horizontal */
					dnorm = dx;
					dpar = dy;
					if((flagField[i][j] & FS_O) == FS_O) {		/* If east neighbor is empty, fluid pressure is from west neighbor */
						Pfluid = P[i-1][j];
					} else if((flagField[i][j] & FS_W) == FS_W) {	/* If west neighbor is empty, fluid pressure is from east neighbor */
						Pfluid = P[i+1][j];
					}
				}
				
				/* Determine distance to neighboring fluid cell */
				m = fmin(fabs(dFdx[i][j]), fabs(dFdy[i][j])) / fmax(fabs(dFdx[i][j]), fabs(dFdy[i][j]));
				Fa = (0.5 * m * dpar * dpar) / dx / dy;
				Fb = 1 - Fa;
				if(fluidFraction[i][j] < Fa) {
					d = 0.5 * dnorm + sqrt(2 * m * fluidFraction[i][j] * dx * dy) - 0.5 * m * dpar;
				} else if(fluidFraction[i][j] < Fb) {
					d = 0.5 * dnorm + fluidFraction[i][j] * dx * dy / dpar;
				} else {
					d = 1.5 * dnorm;				
				}
				eta = dnorm / d;
				
				/* Compute pressure by linear interpolation */
				P[i][j] = (1 - eta) * Pfluid + eta * Psurface;
				P[i][j] = 0;
			}
			
			else if (flagField[i][j] == C_E) {
				P[i][j] = 0;
			}
		}
	}



  /* set boundary values */
  for(i = 1; i <= imax; i++) {
	P[i][0] = P[i][1];
	P[i][jmax+1] = P[i][jmax];
  }
  for(j = 1; j <= jmax; j++) {
	P[0][j] = P[1][j];
	P[imax+1][j] = P[imax][j];
  }  
}

