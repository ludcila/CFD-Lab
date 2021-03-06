#include "boundary_val.h"
#include "helper.h"
#include <string.h>
#include <stdio.h>

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
) {

	int i, j;
	
	/* Left boundary */
	if(il == 1) {
		switch(wl) {
			case BC_NO_SLIP:
				V[0][jb-1] = -V[1][jb-1];
				for(j = jb; j <= jt; j++) {
					U[0][j] = 0;
					V[0][j] = -V[1][j];
				}
				V[0][jt+1] = -V[1][jt+1];
				break;
			case BC_FREE_SLIP:
				V[0][jb-1] = V[1][jb-1];
				for(j = jb; j <= jt; j++) {
					U[0][j] = 0;
					V[0][j] = V[1][j];
				}
				V[0][jt+1] = V[1][jt+1];
				break;
			case BC_OUTFLOW:
				V[0][jb-1] = V[1][jb-1];
				for(j = jb; j <= jt; j++) {
					U[0][j] = U[1][j];
					V[0][j] = V[1][j];
				}
				V[0][jt+1] = V[1][jt+1];
				break;
		}
	}

	/* Right boundary */
	if(ir == imax) {
		switch(wr) {
			case BC_NO_SLIP:
				V[imax+1][jb-1] = -V[imax][jb-1];
				for(j = jb; j <= jt; j++) {
					U[imax][j] = 0;
					V[imax+1][j] = -V[imax][j];
				}
				V[imax+1][jt+1] = -V[imax][jt+1];
				break;
			case BC_FREE_SLIP:
				V[imax+1][jb-1] = V[imax][jb-1];
				for(j = jb; j <= jt; j++) {
					U[imax][j] = 0;
					V[imax+1][j] = V[imax][j];
				}
				V[imax+1][jt+1] = V[imax][jt+1];
				break;
			case BC_OUTFLOW:
				V[imax+1][jb-1] = V[imax][jb-1];
				for(j = jb; j <= jt; j++) {
					U[imax][j] = U[imax-1][j];
					V[imax+1][j] = V[imax][j];
				}
				V[imax+1][jt+1] = V[imax][jt+1];
				break;
		}
	}
	
	/* Top boundary */
	if(jt == jmax) {
		switch(wt) {
			case BC_NO_SLIP:
				U[il-1][jmax+1] = -U[il-1][jmax];
				for(i = il; i <= ir; i++) {
					U[i][jmax+1] = -U[i][jmax];
					V[i][jmax] = 0;
				}
				U[ir+1][jmax+1] = -U[ir+1][jmax];
				break;
			case BC_FREE_SLIP:
				U[il-1][jmax+1] = U[il-1][jmax];
				for(i = il; i <= ir; i++) {
					U[i][jmax+1] = U[i][jmax];
					V[i][jmax] = 0;
				}
				U[ir+1][jmax+1] = U[ir+1][jmax];
				break;
			case BC_OUTFLOW:
				U[il-1][jmax+1] = U[il-1][jmax];
				for(i = il; i <= ir; i++) {
					U[i][jmax+1] = U[i][jmax];
					V[i][jmax] = V[i][jmax-1];
				}
				U[ir+1][jmax+1] = U[ir+1][jmax];
				break;
		}
	}

	/* Bottom boundary */
	if(jb == 1) {
		switch(wb) {
			case BC_NO_SLIP:
				U[il-1][0] = -U[il-1][1];
				for(i = il; i <= ir; i++) {
					U[i][0] = -U[i][1];
					V[i][0] = 0;
				}
				U[ir+1][0] = -U[ir+1][1];
				break;
			case BC_FREE_SLIP:
				U[il-1][0] = U[il-1][1];
				for(i = il; i <= ir; i++) {
					U[i][0] = U[i][1];
					V[i][0] = 0;
				}
				U[ir+1][0] = U[ir+1][1];
				break;
			case BC_OUTFLOW:
				U[il-1][0] = U[il-1][1];
				for(i = il; i <= ir; i++) {
					U[i][0] = U[i][1];
					V[i][0] = V[i][1];
				}
				U[ir+1][0] = U[ir+1][1];
				break;
		}
	}
	
	/* Obstacles */
	for(i = il; i <= ir; i++) {
		for(j = jb; j <= jt; j++) {
		
			switch(Flag[i][j]) {
				case B_N:
					V[i][j] = 0;
					U[i][j] = -U[i][j+1];
					break;
				case B_S:
					V[i][j-1] = 0;
					U[i][j] = -U[i][j-1];
					break;
				case B_O:
					U[i][j] = 0;
					V[i][j] = -V[i+1][j];
					break;
				case B_W:
					U[i-1][j] = 0;
					V[i][j] = -V[i-1][j];
					break;
				case B_NW:
					V[i][j] = 0;
					U[i-1][j] = 0;
					U[i][j] = -U[i][j+1];
					break;
				case B_NO:
					U[i][j] = 0;
					V[i][j] = 0;
					break;
				case B_SW:
					U[i-1][j] = 0;
					V[i][j-1] = 0;
					U[i][j] = -U[i][j-1];
					V[i][j] = -V[i-1][j];
					break;
				case B_SO:
					U[i][j] = 0;
					V[i][j-1] = 0;
					V[i][j] = -V[i+1][j];
					break;
			}
			
		}
	}
	
}

void spec_boundary_val (char *problem, int il, int ir, int jb, int jt, int imax, int jmax, double **U, double **V, double **P, double Re, double xlength, double ylength, double dp){
	
	int i, j;
	
	/* Take care of inflow velocities */
	if(dp == 0) {
		if(strcmp(problem, "flow_over_step") == 0) {
			if(il == 1) {
				for(j = jmax/2 + 1; j <= jt; j++) {
					U[0][j] = 1;
					V[0][j] = 0;
				}
			}
		} else if(strcmp(problem, "karman_vortex_street") == 0) {
			if(il == 1) {
				for(j = jb; j <= jt; j++) {
					U[0][j] = 1;
					V[0][j] = 0;
				}
			}
		} else if(strcmp(problem, "plane_shear_flow") == 0) {
			if(il == 1) {
				for(j = jb; j <= jt; j++) {
					U[0][j] = 1;
					V[0][j] = 0;
				}
			}
		} else if(strcmp(problem, "driven_cavity") == 0) {
			if(jt == jmax) {
				for(i = max(1, il-1); i <= min(imax, ir+1); i++) {
					U[i][jmax + 1] = 2 - U[i][jmax];
					V[i][jmax] = 0;
				}
			}
		} 
	}
	/* Pressure boundary conditions are currently taken care of in the SOR */
	
}
