#include "helper.h"
#include "vof.h"
#include "boundary_val.h"
#include <math.h>

/* Set the initial fluid fraction based on the pgm file */
void init_fluidFraction(int **flagField, double **fluidFraction,double **fluidFraction_alt, int imax, int jmax) {
	int i, j;
	for(i = 1; i <= imax; i++) {
		for(j = 1; j <= jmax; j++) {
			if(flagField[i][j] == C_F) {
				fluidFraction[i][j] = 1;
				fluidFraction_alt[i][j] = 1;
			}
		}
	}
}

/* Bookeeping adjustments as described by Hirt and Nichols */
void adjust_fluidFraction(double **fluidFraction, int **flagField, double epsilon, int imax, int jmax) {

	int i, j;
	int O_empty, W_empty, N_empty, S_empty;
	
	/* TODO: Take care of arbitrary geometry */
	
	for(i = 1; i <= imax; i++) {
		for(j = 1; j <= jmax; j++) {
			
			/* We will consider adjusting and reflagging all cells, except obstacle cells (C_B) */
			if((flagField[i][j] & C_B) == 0) {
			
				O_empty = 0;
				W_empty = 0; 
				N_empty = 0;
				S_empty = 0;
			
				/* Treat empty cells */
				if(fluidFraction[i][j] < epsilon) {
					fluidFraction[i][j] = 0;
					flagField[i][j] = C_E;
				}
			
				/* Treat cells that contain fluid (either full or a fraction) */
				else {
			
					/* Initially set as full fluid cell (no empty neighbors) */
					flagField[i][j] = C_F;
				
					/* Determine empty neighbors with three conditions:
						(neighbor is not domain boundary) AND (is not obstacle) AND (has less than eps fluid) */
					O_empty = (i+1 <= imax)	&& !(flagField[i+1][j] & C_B) 	&& (fluidFraction[i+1][j] < epsilon);
					W_empty = (i-1 >= 1)	&& !(flagField[i-1][j] & C_B)	&& (fluidFraction[i-1][j] < epsilon);
					N_empty = (j+1 <= jmax)	&& !(flagField[i][j+1] & C_B)	&& (fluidFraction[i][j+1] < epsilon);
					S_empty = (j-1 >= 1)	&& !(flagField[i][j-1] & C_B)	&& (fluidFraction[i][j-1] < epsilon);
				
					/* Toggle bit for each empty neighbor empty neighbors
						if there are no empty neighbors, the flag will remain unchanged (i.e., C_F)
						if there are empty neighbors, the flag will be a combination of FS_x flags
					*/
					flagField[i][j] = flagField[i][j] | (O_empty * FS_O) | (W_empty * FS_W) | (N_empty * FS_N) | (S_empty * FS_S);
				
					/* If it is a full fluid cell, we reduce the fluid fraction by 1.1*epsilon for each empty neighbor */
					if(fluidFraction[i][j] > 1-epsilon) {
						fluidFraction[i][j] = 1 - 1.1 * epsilon * (O_empty + W_empty + N_empty + S_empty);
					}
					
				}
				
			}
			
		}
	}

}

/* Determine orientation of the free surface */
void calculate_freeSurfaceOrientation(double **fluidFraction, int **flagField, double **dFdx, double **dFdy, double dx, double dy, int imax, int jmax) {

	int i, j;
	
	/* Copy fluidFraction to boundary cells, just to make computations easier */
	for(i = 1; i <= imax; i++) {
		fluidFraction[i][0] = fluidFraction[i][1];
		fluidFraction[i][jmax+1] = fluidFraction[i][jmax];
	}
	for(j = 1; j <= jmax; j++) {
		fluidFraction[0][j] = fluidFraction[1][j];
		fluidFraction[imax+1][j] = fluidFraction[imax][j];
	}
	
	/* (Assuming uniform grid in each axis) */
	for(i = 1; i <= imax; i++) {
		for(j = 1; j <= jmax; j++) {
			if(flagField[i][j] & C_FS) {
				dFdx[i][j] = 0.5 * ((fluidFraction[i+1][j-1] + fluidFraction[i+1][j] + fluidFraction[i+1][j+1]) - (fluidFraction[i-1][j-1] + fluidFraction[i-1][j] + fluidFraction[i-1][j+1])) / dx;
				dFdy[i][j] = 0.5 * ((fluidFraction[i-1][j+1] + fluidFraction[i][j+1] + fluidFraction[i+1][j+1]) - (fluidFraction[i-1][j-1] + fluidFraction[i][j-1] + fluidFraction[i+1][j-1])) / dy;
			}
		}
	}
	
}

/* Timestepping for the fluid fraction field */
void calculate_fluidFraction(
	double **fluidFraction,
 	double **fluidFraction_alt,
 	int **flagField,
	double **U, 
	double **V, 
	double **dFdx, 
	double **dFdy, 
	int imax,
	int jmax,
	double dx, 
	double dy, 
	double dt
) {

	int i,j,sign;
	double F_D_x,F_AD_x;
	double V_x,CF_x,delta_F_right,delta_F_left;
	double F_D_y,F_AD_y;
	double V_y,CF_y,delta_F_bottom,delta_F_top;

  
	for(i=1;i<=imax;i++){
		for(j=1;j<=jmax;j++){
			fluidFraction_alt[i][j]=fluidFraction[i][j];
		}
	}

	/* ======================= X sweep ====================== */
	
    for (j=1;j<=jmax;j++){
        for (i=1;i<=imax-1;i++){
		        	
			if(i==1){
				delta_F_left=0;
			}
		
			/* Right going */
			if(U[i][j]>0){
				if((i < imax && fluidFraction_alt[i+1][j] == 0) || (i > 1 && fluidFraction_alt[i-1][j] == 0)) {
					F_AD_x = fluidFraction_alt[i+1][j];
				} else if((flagField[i][j] & C_FS) == C_FS){
					if(fabs(dFdy[i][j]) > fabs(dFdx[i][j])) {
						F_AD_x = fluidFraction_alt[i][j];
					} else {
						F_AD_x = fluidFraction_alt[i+1][j];
					}
				} else {
					F_AD_x = fluidFraction_alt[i][j];
				}
				F_D_x=fluidFraction_alt[i][j];
				sign=-1;
			}
		
			/* Left going */
			else if(U[i][j]<0){
				if(fluidFraction_alt[i][j] == 0 || (i < imax-1 && fluidFraction_alt[i+2][j] == 0)) {
					F_AD_x = fluidFraction_alt[i][j];
				} else if((flagField[i][j] & C_FS) == C_FS){
					if(fabs(dFdy[i][j]) > fabs(dFdx[i][j])) {
						F_AD_x = fluidFraction_alt[i+1][j];
					} else {
						F_AD_x = fluidFraction_alt[i][j];
					}
				} else {
					F_AD_x = fluidFraction_alt[i+1][j];
				}
				F_D_x=fluidFraction_alt[i+1][j];
				sign=1;
			}
		
			/* Zero velocity */
			else {
				F_D_x=0;
			}

			/* Update fluid fraction based on computed fluxes */
		    V_x=U[i][j]*dt;                
		    CF_x=fmax((1.0-F_AD_x)*fabs(V_x)-(1.0-F_D_x)*dx,0.0);
		    delta_F_right=fmin(F_AD_x*fabs(V_x)+CF_x,F_D_x*dx)/dx;
			fluidFraction[i][j]=fluidFraction[i][j]+sign*delta_F_right+delta_F_left;
			delta_F_left=(-1)*sign*delta_F_right;
		
		 }
		
		/* Treat right-most cell */
		fluidFraction[imax][j]=fluidFraction[imax][j]+delta_F_left;
	
	}
  
	/* ======================= Swap intermediate solution ====================== */
	for(i=1;i<=imax;i++){
		for(j=1;j<=jmax;j++){
			fluidFraction_alt[i][j]=fluidFraction[i][j];
		}
	}
 
	/* ======================= Y sweep ====================== */
    for (i=1;i<=imax;i++){
        for (j=1;j<=jmax-1;j++){

			if(j==1){
				delta_F_bottom=0;
			}

			/* Up going */
			if(V[i][j]>0){
				if((j < jmax && fluidFraction_alt[i][j+1] == 0) || (j > 1 && fluidFraction_alt[i][j-1] == 0)) {
					F_AD_y = fluidFraction_alt[i][j+1];
				}else if((flagField[i][j] & C_FS) == C_FS){
					if(fabs(dFdy[i][j]) < fabs(dFdx[i][j])) {
						F_AD_y = fluidFraction_alt[i][j];
					} else {
						F_AD_y = fluidFraction_alt[i][j+1];
					}
				} else {
					F_AD_y = fluidFraction_alt[i][j];
				}
				F_D_y=fluidFraction_alt[i][j];
				sign=-1;
			
			}

			/* Down going */
			else if(V[i][j]<0){
				if(fluidFraction_alt[i][j] == 0 || (j < jmax-1 && fluidFraction_alt[i][j+2] == 0)) {
					F_AD_y = fluidFraction_alt[i][j];
				}else if((flagField[i][j] & C_FS) == C_FS){
					if(fabs(dFdy[i][j]) < fabs(dFdx[i][j])) {
						F_AD_y = fluidFraction_alt[i][j+1];
					} else {
						F_AD_y = fluidFraction_alt[i][j];
					}
				} else {
					F_AD_y = fluidFraction_alt[i][j+1];
				}
				F_D_y=fluidFraction_alt[i][j+1];
				sign=1;
			}
		
			/* Zero velocity */
			else {
				F_D_y=0;
			}
		
			/* Update fluid fraction based on computed fluxes */
		    V_y=V[i][j]*dt;                
		    CF_y=fmax((1.0-F_AD_y)*fabs(V_y)-(1.0-F_D_y)*dy,0.0);
		    delta_F_top=fmin(F_AD_y*fabs(V_y)+CF_y,F_D_y*dy)/dy;
			fluidFraction[i][j]=fluidFraction[i][j]+sign*delta_F_top+delta_F_bottom;
			delta_F_bottom=(-1)*sign*delta_F_top;
				
		}

		/* Treat top-most cell */
		fluidFraction[i][jmax]=fluidFraction[i][jmax]+delta_F_bottom;
	
	}


}
