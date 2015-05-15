#include "streaming.h"
#include "LBDefinitions.h"
#include <stdlib.h>

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
  	int x,y,z;
	int i;

	for(z=1;z<=xlength;z++){
		for(y=1;y<=xlength;y++){
			for(x=1;x<=xlength;x++){
/*				if(flagField[(z*xlength*xlength+y*xlength+x)]==0){*/
					for(i=0;i<Q;i++){
						streamField[Q*(z*(xlength+1)*(xlength+1)+y*(xlength+1)+x)+i]=collideField
							[Q*((-LATTICEVELOCITIES[i][1]+z)*(xlength+1)*(xlength+1)+
							(-LATTICEVELOCITIES[i][2]+y)*(xlength+1)+(-LATTICEVELOCITIES[i][3]+x))+i]; /*need to be checked*/
							}
			/*	}    */
			}	
		}
	}
}

