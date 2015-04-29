#include "boundary_val.h"

void boundaryvalues(
	int imax,
	int jmax,
	double **U,
	double **V
	
) {
	/* Modified by San Yu *//*assume grid from 0->imax, so U in (imax+1)*(jmax+2)  */
	int i,j;	    /* V in (imax+2)*(jmax+1) */	    	
	double U_wall=10;   /*needs to be initialized*/

	for(j=1;j<=jmax;j++){	/*left wall*/
		U[0][j]=0;
		V[0][j]=-V[1][j];
	}
	for(j=1;j<=jmax;j++){  /*right wall*/
		U[imax][j]=0;
		V[imax+1][j]=-V[imax][j];
	}
	for(i=1;i<=imax;i++){	/*floor*/
		U[i][0]=-U[i][1];
		V[i][0]=0;
	}
	for(i=1;i<=imax;i++){	/*ceiling*/
		U[i][jmax+1]=-U[i][jmax]+2*U_wall;
		V[i][jmax]=0;
	}

}
