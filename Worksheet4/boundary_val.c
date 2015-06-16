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
	
	/* Extended indices (e.g. compute U from il-1 to ir+1, unless il=1 or ir=imax) */
	int il_ext = max(1, il-1);
	int ir_ext = min(imax, ir+1);
	int jb_ext = max(1, jb-1);
	int jt_ext = min(jmax, jt+1);
	
	/* Left boundary */
	if(il == 1) {
		switch(wl) {
			case BC_NO_SLIP:
				for(j = jb; j <= jt; j++) {
					U[0][j] = 0;
					V[0][j] = -V[1][j];
				}
				/* Extended range for V */
				V[0][jb_ext] = -V[1][jb_ext];
				V[0][jt_ext] = -V[1][jt_ext];
				break;
			case BC_FREE_SLIP:
				for(j = jb; j <= jt; j++) {
					U[0][j] = 0;
					V[0][j] = V[1][j];
				}
				break;
			case BC_OUTFLOW:
				for(j = jb; j <= jt; j++) {
					U[0][j] = U[1][j];
					V[0][j] = V[1][j];
				}
				break;
		}
	}

	/* Right boundary */
	if(ir == imax) {
		switch(wr) {
			case BC_NO_SLIP:
				for(j = jb; j <= jt; j++) {
					U[imax][j] = 0;
					V[imax+1][j] = -V[imax][j];
				}
				/* Extended range for V */
				V[imax+1][jb_ext] = -V[imax][jb_ext];
				V[imax+1][jt_ext] = -V[imax][jt_ext];
				break;
			case BC_FREE_SLIP:
				for(j = jb; j <= jt; j++) {
					U[imax][j] = 0;
					V[imax+1][j] = V[imax][j];
				}
				break;
			case BC_OUTFLOW:
				for(j = jb; j <= jt; j++) {
					U[imax][j] = U[imax-1][j];
					V[imax+1][j] = V[imax][j];
				}
				break;
		}
	}
	
	/* Top boundary */
	if(jt == jmax) {
		switch(wt) {
			case BC_NO_SLIP:
				for(i = il; i <= ir; i++) {
					U[i][jmax+1] = -U[i][jmax];
					V[i][jmax] = 0;
				}
				/* Extended range for U */
				U[il_ext][jmax+1] = -U[il_ext][jmax];
				U[ir_ext][jmax+1] = -U[ir_ext][jmax];
				break;
			case BC_FREE_SLIP:
				for(i = il; i <= ir; i++) {
					U[i][jmax+1] = U[i][jmax];
					V[i][jmax] = 0;
				}
				break;
			case BC_OUTFLOW:
				for(i = il; i <= ir; i++) {
					U[i][jmax+1] = U[i][jmax];
					V[i][jmax] = V[i][jmax-1];
				}
				break;
		}
	}

	/* Bottom boundary */
	if(jb == 1) {
		switch(wb) {
			case BC_NO_SLIP:
				for(i = il; i <= ir; i++) {
					U[i][0] = -U[i][1];
					V[i][0] = 0;
				}
				/* Extended range for U */
				U[il_ext][0] = -U[il_ext][1];
				U[ir_ext][0] = -U[ir_ext][1];
				break;
			case BC_FREE_SLIP:
				for(i = il; i <= ir; i++) {
					U[i][0] = U[i][1];
					V[i][0] = 0;
				}
				break;
			case BC_OUTFLOW:
				for(i = il; i <= ir; i++) {
					U[i][0] = U[i][1];
					V[i][0] = V[i][1];
				}
				break;
		}
	}
	
	/* Obstacles */
	for(i = il_ext; i <= ir_ext; i++) {
		for(j = jb_ext; j <= jt_ext; j++) {
		
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
