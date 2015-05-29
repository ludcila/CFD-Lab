#include "boundary_val.h"

/* flags
* C_F: fluid cell
* C_B: obstacle cell
* B_xy: boundary cell
*/

enum B_xy {B_NULL = 0, B_N = 1, B_S = 2, B_W = 4, B_NW = 5, B_SW = 6, B_NSW = 7, B_O = 8, B_NO = 9, B_SO = 10, B_NSO = 11, B_NWO = 13, B_SWO = 14};


/*
* wl, wr, wt, wb 
* 1 = noslip
* 2 = feeslip
* 3 = outflow
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
) {

	int i, j;	
	
	if (wl == 1){
		noslip(imax, jmax, U, V, 0);
	}
	else if(wl == 2){
		freeslip(imax, jmax, U, V, 0);
	}
	else if(wl == 3){
		outflow(imax, jmax, U, V, 0);
	}

	if (wr == 1){
		noslip(imax, jmax, U, V, 1);
	}
	else if(wr == 2){
		freeslip(imax, jmax, U, V, 1);
	}
	else if(wr == 3){
		outflow(imax, jmax, U, V, 1);
	}

	if (wt == 1){
		noslip(imax, jmax, U, V, 2);
	}
	else if(wt == 2){
		freeslip(imax, jmax, U, V, 2);
	}
	else if(wt == 3){
		outflow(imax, jmax, U, V, 2);
	}


	if (wb == 1){
		noslip(imax, jmax, U, V, 3);
	}
	else if(wb == 2){
		freeslip(imax, jmax, U, V, 3);
	}
	else if(wb == 3){
		outflow(imax, jmax, U, V, 3);
	}


	/* set obstacles*/
	/* treat v */
	for (i = 1; i <= imax - 1; i++) {
		for (j = 1; j <= jmax; j++) {
			if(Flag[i][j] && 1){/*B_N, north cell is fluid*/
				V[i][j] = 0;
			} else if(!(Flag[i][j] && 1)){
				V[i][j] = -V[i-1][j];
			}

			if(Flag[i][j] && 2){/*B_S, south cell is fluid*/
				V[i][j-1] = 0;	
			} else if(!(Flag[i][j] && 2)){
				V[i][j-1] = -V[i-1][j-1];
			}

		}
	}
	/* treat u */
	for (i = 1; i <= imax; i++) {
		for (j = 1; j <= jmax - 1; j++) {
			if(Flag[i][j] && 4){/*B_W, west cell is fluid*/
				U[i-1][j] = 0;
			} else if(!(Flag[i][j] && 4)){
				U[i-1][j] = -U[i-1][j+1];	
			}

			if(Flag[i][j] && 8){/*B_O, east cell is fluid*/
				U[i][j] = 0;
			} else if(!(Flag[i][j] && 8)){
				U[i][j] = -U[i][j+1];
			}
			
		}
	}



}

/**
* side = 0 -> left wall, side = 1 -> right wall, side = 2 -> top wall, side = 3 -> bottom wall
*
*/
void noslip(int imax, int jmax, double **U, double **V, int side){
	
	int i, j;

	if(side == 0){
		/* Left wall */	
		for(j = 0; j <= jmax; j++) {
			U[0][j] = 0;
			V[0][j] = -V[1][j];
		}
	}else if(side == 1){
		/* Right wall */
		for(j = 0; j <= jmax; j++) {
			U[imax][j] = 0;
			V[imax+1][j] = -V[imax][j];
		}
	} else if (side == 2){		
		/* Top */
		for(i = 0; i <= imax; i++) {
			U[i][jmax+1] = -U[i][jmax];
			V[imax][0] = 0;
		}
	} else if (side == 3){
		/* Bottom */
		for(i = 0; i <= imax; i++) {
			U[i][0] = -U[i][1];
			V[i][0] = 0;
		}
	}
}

void freeslip(int imax, int jmax, double **U, double **V, int side){
	
	int i, j;

	if(side == 0){
		/* Left wall */	
		for(j = 0; j <= jmax; j++) {
			U[0][j] = 0;
			V[0][j] = V[1][j];
		}
	}else if(side == 1){
		/* Right wall */
		for(j = 0; j <= jmax; j++) {
			U[imax][j] = 0;
			V[imax+1][j] = V[imax][j];
		}
	} else if (side == 2){		
		/* Top */
		for(i = 0; i <= imax; i++) {
			U[i][jmax+1] = U[i][jmax];
			V[imax][0] = 0;
		}
	} else if (side == 3){
		/* Bottom */
		for(i = 0; i <= imax; i++) {
			U[i][0] = U[i][1];
			V[i][0] = 0;
		}
	}
}

void outflow(int imax, int jmax, double **U, double **V, int side){
	
	int i, j;

	if(side == 0){
		/* Left wall */	
		for(j = 0; j <= jmax; j++) {
			U[0][j] = U[1][j];
			V[0][j] = V[1][j];
		}
	}else if(side == 1){
		/* Right wall */
		for(j = 0; j <= jmax; j++) {
			U[imax][j] = V[imax-1][j];
			V[imax+1][j] = V[imax][j];
		}
	} else if (side == 2){		
		/* Top */
		for(i = 0; i <= imax; i++) {
			U[i][jmax+1] = U[i][jmax];
			V[imax][0] = V[imax][1];
		}
	} else if (side == 3){
		/* Bottom */
		for(i = 0; i <= imax; i++) {
			U[i][0] = U[i][1];
			V[i][0] = V[i][1];
		}
	}
}

void movingwall(int imax, int jmax, double **U, double **V, int side){

	int i, j;
	double U_wall = 1.0;
	if (side == 0){
		/* Left wall */
		for(j = 0; j <= jmax; j++) {
			U[0][j] = 0;
			V[0][j] = 2 * U_wall - V[1][j];
		}
	} else if (side == 1){
		/* Right wall */
		for(j = 0; j <= jmax; j++) {
			U[imax][j] = 0;
			V[imax+1][j] = 2 * U_wall - V[imax][j];
		}
	} else if (side == 2){
		/* Top */
		for(i = 0; i <= imax; i++) {
			U[i][jmax+1] = 2 * U_wall - U[i][jmax];
			V[i][jmax] = 0;
		}
	} else if (side == 3){
		/* Bottom */
		for(i = 0; i <= imax; i++) {
			U[i][0] = 2 * U_wall - U[i][1];
			V[i][0] = 0;
		}
	}
}

void spec_boundary_val (char *problem, int imax, int jmax, double **U, double **V){

}
