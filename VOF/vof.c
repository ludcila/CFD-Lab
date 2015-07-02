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
				O_empty = fluidFraction[min(i+1, imax)][j] < epsilon;
				W_empty = fluidFraction[max(i-1, 1)][j] < epsilon;
				N_empty = fluidFraction[i][min(j+1, jmax)] < epsilon;
				S_empty = fluidFraction[i][max(j-1, 1)] < epsilon;
				
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
	double F_D_right,F_D_left,F_AD_left,F_AD_right;
	double V_x_right,V_x_left,CF_x_right,CF_x_left,delta_F_right,delta_F_left;
	double F_D_bottom,F_D_top,F_AD_top,F_AD_bottom,V_y_bottom,V_y_top,CF_y_bottom,CF_y_top;
	double delta_F_bottom,delta_F_top;




/*=======================for x sweep======================*/
    for (i=1;i<=imax;i++){
        for (j=1;j<=jmax;j++){
            
                if (j<jmax){
                    if (U[i][j]>0){
                        F_D_right=fluidFraction_alt[i][j];
                        sign=1;}
                    else if(U[i][j]<0){
                        F_D_right=fluidFraction_alt[i][j+1];
                        sign=-1;
                    }
                }
                if (j>1){
                    if(U[i][j-1]>0){
                        F_D_left=fluidFraction_alt[i][j-1];
                        sign=1;}
                    else if(U[i][j-1]<0){
                        F_D_left=fluidFraction_alt[i][j];
                        sign=-1;
                    }
                }
                
                if(j==1){
                    if(U[i][j]>0)
                        F_D_left=0.0;
                    else
                        F_D_left=fluidFraction_alt[i][j];
                    }
                else if(j==jmax){
                    if(U[i][j]<0)
                        F_D_right=0.0;
                    else
                        F_D_right=fluidFraction_alt[i][j];
                    
                }

                 F_AD_left=F_D_left;
                 F_AD_right=F_D_right;

         
                
                V_x_right=U[i][j]*dt;
                V_x_left=U[i-1][j]*dt;
                CF_x_right=fmax((1.0-F_AD_right)*abs(V_x_right)-(1.0-F_D_right)*dx,0);
                CF_x_left=fmax((1.0-F_AD_left)*abs(V_x_left)-(1.0-F_D_left)*dx,0);
                delta_F_right=fmin(F_AD_right*abs(V_x_right)+CF_x_right,F_D_right*dx)*dx;
              	delta_F_left=fmin(F_AD_left*abs(V_x_left)+CF_x_left,F_AD_left*dx)*dx;
                fluidFraction[i][j]=fluidFraction[i][j]-sign*delta_F_right+sign*delta_F_left;
            
    
         }
    }
  

   
/*=======================for y sweep======================*/
    for (j=1;j<=jmax;j++){
        for (i=1;i<=imax;i++){

               if (i<imax){
                    if (V[i][j]>0){
                        F_D_bottom=fluidFraction_alt[i][j];
                        sign=1;}
                    else if(V[i][j]<0){
                        F_D_bottom=fluidFraction_alt[i+1][j];
                        sign=-1;
                    }
                }
                if (i>1){
                    if((V[i-1][j]>0)){
                        F_D_top=fluidFraction_alt[i-1][j];
                        sign=1;}
                    else if(V[i-1][j]<0){
                        F_D_top=fluidFraction_alt[i][j];
                        sign=-1;
                    }
                }
                
                if(i==1){
                    if(V[i][j]>0)
                        F_D_top=0.0;
                    else
                        F_D_top=fluidFraction_alt[i][j];
                    }
                else if(i==imax){
                    if(V[i][j]<0)
                        F_D_bottom=0.0;
                    else
                        F_D_bottom=fluidFraction_alt[i][j];
                   
                }

                 F_AD_top=F_D_top;
                 F_AD_bottom=F_D_bottom;


                V_y_bottom=V[i][j-1]*dt;
                V_y_top=V[i][j]*dt;
                CF_y_bottom=fmax((1.0-F_AD_bottom)*abs(V_y_bottom)-(1.0-F_D_bottom)*dy,0);
                CF_y_top=fmax((1.0-F_AD_top)*abs(V_y_top)-(1.0-F_D_top)*dy,0);
                delta_F_bottom=fmin(F_AD_bottom*abs(V_y_bottom)+CF_y_bottom,F_D_bottom*dy)*dy;
                delta_F_top=fmin(F_AD_top*abs(V_y_top)+CF_y_top,F_AD_top*dy)*dy;
                fluidFraction[i][j]=fluidFraction[i][j]-sign*delta_F_bottom+sign*delta_F_top;
    
         }
    }


    
	for(i=1;i<=imax;i++){
		for(j=1;j<=jmax;j++){
			fluidFraction_alt[i][j]=fluidFraction[i][j];
			}
	}










}
