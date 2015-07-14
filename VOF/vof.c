#include "helper.h"
#include "vof.h"
#include <math.h>

/* Set the initial fluid fraction based on the pgm file */
void init_fluidFraction(int **pgm, double **fluidFraction,double **fluidFraction_alt, int imax, int jmax) {
	int i, j;
	for(i = 1; i <= imax; i++) {
		for(j = 1; j <= jmax; j++) {
			if(pgm[i][j] == 1) {
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
				
				/* Determine empty neighbors; the min/max prevent it from checking neighbors on the domain boundary */
				if(i+1 <= imax) {
					O_empty = fluidFraction[i+1][j] < epsilon;
				}
				if(i-1 >= 1) {
					W_empty = fluidFraction[i-1][j] < epsilon;
				}
				if(j+1 <= jmax) {
					N_empty = fluidFraction[i][j+1] < epsilon;
				}
				if(j-1 >= 1) {
					S_empty = fluidFraction[i][j-1] < epsilon;
				}
				
				/* Toggle bit for each empty neighbor empty neighbors
					if there are no empty neighbors, the flag will remain unchanged (i.e., C_F)
					if there are empty neighbors, the flag will be a combination of FS_x flags
				*/
				flagField[i][j] = flagField[i][j] | (O_empty * FS_O) | (W_empty * FS_W) | (N_empty * FS_N) | (S_empty * FS_S);
				
				/* If it is a full fluid cell, we reduce the fluid fraction by 1.1*epsilon for each empty neighbor */
				if(fluidFraction[i][j] > 1-epsilon)
					fluidFraction[i][j] = 1 - 1.1 * epsilon * (O_empty + W_empty + N_empty + S_empty);
				
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


/*=======================for x sweep======================*/
    for (j=1;j<=jmax;j++){
        for (i=1;i<=imax-1;i++){
            	
		if(i==1){
			delta_F_left=0;
		}

		if(U[i][j]>0){
/*<<<<<<< HEAD
			if(fluidFraction_alt[i+1][j]==0 || fluidFraction_alt[i-1][j]==0){
=======*/
			if(fluidFraction_alt[i+1][j] ==0 || fluidFraction_alt[i-1][j] < 1e-6){
/*>>>>>>> 5df4228685e99ead1c7676d3ca1b9fa4acbe9fb4*/
				F_AD_x = fluidFraction_alt[i+1][j];
			}else{
				if(fabs(dFdy[i][j]) > fabs(dFdx[i][j])) {
					F_AD_x = fluidFraction_alt[i][j];
				} else {
					F_AD_x = fluidFraction_alt[i+1][j];
				}
			}
			F_D_x=fluidFraction_alt[i][j];
			sign=-1;
		}else if(U[i][j]<0){
			if(fluidFraction_alt[i][j] < 1e-6 || fluidFraction_alt[i+2][j] < 1e-6){
				F_AD_x = fluidFraction_alt[i][j];
			}else{
				if(fabs(dFdy[i][j]) > fabs(dFdx[i][j])) {
					F_AD_x = fluidFraction_alt[i+1][j];
				} else {
					F_AD_x = fluidFraction_alt[i][j];
				}
			}
			F_D_x=fluidFraction_alt[i+1][j];
			sign=1;}
		else
			F_D_x=0;

	
			         

                V_x=U[i][j]*dt/2;                
                CF_x=fmax((1.0-F_AD_x)*fabs(V_x)-(1.0-F_D_x)*dx,0.0);
/*<<<<<<< HEAD
                delta_F_right=fmin(F_AD_x*fabs(V_x)+CF_x,F_D_x*dx);
                fluidFraction[i][j]=fluidFraction[i][j]+sign*delta_F_right+delta_F_left;
=======*/
                delta_F_right=fmin(F_AD_x*fabs(V_x)+CF_x,F_D_x*dx)/dx;
                /*fluidFraction[i][j]=fluidFraction_alt[i][j]+sign*delta_F_right+delta_F_left;*/
		fluidFraction[i][j]=fluidFraction[i][j]+sign*delta_F_right+delta_F_left;
/*>>>>>>> 5df4228685e99ead1c7676d3ca1b9fa4acbe9fb4*/

    		delta_F_left=(-1)*sign*delta_F_right;
         }


	/*fluidFraction[imax][j]=fluidFraction_alt[imax][j]+delta_F_left;*/
	fluidFraction[imax][j]=fluidFraction[imax][j]+delta_F_left;
    }
  
	for(i=1;i<=imax;i++){
		for(j=1;j<=jmax;j++){
			fluidFraction_alt[i][j]=fluidFraction[i][j];
			}
	}
 
/*=======================for y sweep======================*/
    for (i=1;i<=imax;i++){
        for (j=1;j<=jmax-1;j++){

               if(j==1){
			delta_F_bottom=0;
		}


		if(V[i][j]>0){
/*<<<<<<< HEAD
			if(fluidFraction_alt[i][j+1]==0 || fluidFraction_alt[i][j-1]==0){
=======*/
			if(fluidFraction_alt[i][j+1] < 1e-6 || fluidFraction_alt[i][j-1] < 1e-6){
/*>>>>>>> 5df4228685e99ead1c7676d3ca1b9fa4acbe9fb4*/
				F_AD_y = fluidFraction_alt[i][j+1];
			}else{
				if(fabs(dFdy[i][j]) < fabs(dFdx[i][j])) {
					F_AD_y = fluidFraction_alt[i][j];
				} else {
					F_AD_y = fluidFraction_alt[i][j+1];
				}
			}
			F_D_y=fluidFraction_alt[i][j];
			sign=-1;}
		else if(V[i][j]<0){
					
			if(fluidFraction_alt[i][j] < 1e-6 || fluidFraction_alt[i][j+2] < 1e-6){
				F_AD_y = fluidFraction_alt[i][j];
			}else{	
				if(fabs(dFdy[i][j]) < fabs(dFdx[i][j])) {
					F_AD_y = fluidFraction_alt[i][j+1];
				} else {
					F_AD_y = fluidFraction_alt[i][j];
				}
			}
			F_D_y=fluidFraction_alt[i][j+1];
			sign=1;}
		else
			F_D_y=0;
       
/*                F_AD_y=F_D_y;*/
                V_y=V[i][j]*dt/2;                
                CF_y=fmax((1.0-F_AD_y)*fabs(V_y)-(1.0-F_D_y)*dy,0.0);
/*<<<<<<< HEAD
                delta_F_top=fmin(F_AD_y*fabs(V_y)+CF_y,F_D_y*dy);
                fluidFraction[i][j]=fluidFraction[i][j]+sign*delta_F_top+delta_F_bottom;
=======*/
                delta_F_top=fmin(F_AD_y*fabs(V_y)+CF_y,F_D_y*dy)/dy;
                /*fluidFraction[i][j]=fluidFraction_alt[i][j]+sign*delta_F_top+delta_F_bottom;*/
fluidFraction[i][j]=fluidFraction[i][j]+sign*delta_F_top+delta_F_bottom;
/*>>>>>>> 5df4228685e99ead1c7676d3ca1b9fa4acbe9fb4*/
            	
    		delta_F_bottom=(-1)*sign*delta_F_top;
    		
         }

	/*fluidFraction[i][jmax]=fluidFraction_alt[i][jmax]+delta_F_bottom;*/
fluidFraction[i][jmax]=fluidFraction[i][jmax]+delta_F_bottom;
    }


    
	for(i=1;i<=imax;i++){
		for(j=1;j<=jmax;j++){
			fluidFraction_alt[i][j]=fluidFraction[i][j];
			}
	}










}
