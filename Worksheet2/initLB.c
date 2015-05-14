#include "initLB.h"
#include <math.h>
#include "LBDefinitions.h"

int readParameters( const char *szFileName, int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){


	READ_INT(szFileName,*xlength);
	READ_DOUBLE(szFileName,*tau);
	READ_DOUBLE(szFileName,*velocityWall);
	READ_INT(szFileName,*timesteps);
	READ_INT(szFileName,*timestepsPerPlotting);
	


  return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
 

	int x,y,z,i;

/* Initialize flageField: geometry mapping for a cell: FLUID=0, NO SLIP=1 and MOVING WALL=2 */	
	for(z=0;z<=xlength+1;z++){
		for(y=0;y<=xlength+1;y++){
			for(x=0;x<=xlength+1;x++){
				if(x==0 || z==0 || z==xlength+1){
					flagField[z*(xlength+1)*(xlength+1)+y*(xlength+1)+x]=1;}
				else if(x==xlength+1){
					flagField[z*(xlength+1)*(xlength+1)+y*(xlength+1)+x]=2;}
				else if(y==0 || y==xlength+1){
					flagField[z*(xlength+1)*(xlength+1)+y*(xlength+1)+x]=1;}
				else{
					flagField[z*(xlength+1)*(xlength+1)+y*(xlength+1)+x]=0;}
			}
		}
	}
/* Initialize matrice of collision and stream:f(x,t=0)=omega_i*/
	for(z=1;z<=xlength;z++){
		for(y=1;y<=xlength;y++){
			for(x=1;x<=xlength;x++){
				for(i=0;i<Q;i++){
					if(i==9){
						collideField[Q*(z*(xlength+1)*(xlength+1)+y*(xlength+1)+x)+i]=W_12_36;
						streamField[Q*(z*(xlength+1)*(xlength+1)+y*(xlength+1)+x)+i]=W_12_36;}
					else if(i==2 ||i==6 ||i==8 ||i==10 ||i==12 ||i==16){
						collideField[Q*(z*(xlength+1)*(xlength+1)+y*(xlength+1)+x)+i]=W_2_36;
						streamField[Q*(z*(xlength+1)*(xlength+1)+y*(xlength+1)+x)+i]=W_2_36;}
					else{
						collideField[Q*(z*(xlength+1)*(xlength+1)+y*(xlength+1)+x)+i]=W_1_36;
						streamField[Q*(z*(xlength+1)*(xlength+1)+y*(xlength+1)+x)+i]=W_1_36;}
						}
					      }
					}
				}


}

