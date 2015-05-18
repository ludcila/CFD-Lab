#include "helper.h"
#include "visualLB.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdio.h>

void writeVtkOutput(const double * const collideField, const int * const flagField, const char * filename, unsigned int tau, int xlength) {
	/* TODO Jae*/
	int i, x, y, z;
	double density;
	char szFileName[80];
	int cellIdx;
	int numGridPoints = xlength + 2;
	FILE *fp=NULL;

	sprintf( szFileName, "%s_%s.vtk", filename, "CFD_group1"  );

	if( fp == NULL )		       
	{
		char szBuff[80];
		sprintf( szBuff, "Failed to open %s", szFileName );
		ERROR( szBuff );
		return;
	}
	write_vtkHeader( fp, xlength );
	write_vtkPointCoordinates(fp, xlength, tau);

	fprintf(fp,"POINT_DATA %i \n", (xlength+1)*(xlength+1) );

	fprintf(fp,"\n");
	fprintf(fp, "density\n");
	
	for(z = 1; z <= xlength; z++) {
		for(y = 1; y <= xlength; y++) {
			for(x = 1; x <= xlength; x++) {
			
				cellIdx = Q * (z * numGridPoints * numGridPoints + y * numGridPoints + x);

/*
				const double *const currentCell = collideField + cellIdx;

				for(i = 0; i < Q; i++) {
					density += currentCell[i];
				}
*/
				computeDensity(collideField + cellIdx, &density);
				fprintf(fp, "%f\n", density);
			}
		}
	}

	fprintf(fp,"\n");
	fprintf(fp, "velocity (X, Y, Z)\n");
	for(i = 0; i < Q; i++) {
		fprintf(fp, "%d %d %d\n", LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]);
	}

	if( fclose(fp) )
	{
		char szBuff[80];
		sprintf( szBuff, "Failed to close %s", szFileName );
		ERROR( szBuff );
	}
}

void write_vtkHeader( FILE *fp, int xlength ) {
	if( fp == NULL )		       
	{
		char szBuff[80];
		sprintf( szBuff, "Null pointer in write_vtkHeader" );
		ERROR( szBuff );
		return;
	}

	fprintf(fp,"# vtk DataFile Version 2.2\n");
	fprintf(fp,"generated by CFD-lab course output (written by Tobias Neckel edited by Group 1) \n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"\n");	
	fprintf(fp,"DATASET STRUCTURED_GRID\n");
	fprintf(fp,"DIMENSIONS  %i %i %i \n", xlength+1, xlength+1, xlength+1);
	fprintf(fp,"POINTS %i float\n", (xlength+1)*(xlength+1)*(xlength+1) );
	fprintf(fp,"\n");
}


void write_vtkPointCoordinates( FILE *fp, int xlength, double tau) {
	double originX = 0.0;  
	double originY = 0.0;
	double originZ = 0.0;

	int i = 0;
	int j = 0;
	int k = 0;

	for(j = 0; j < xlength+1; j++) {
		for(i = 0; i < xlength+1; i++) {
			fprintf(fp, "%f %f %f\n", originX+(i*tau), originY+(j*tau), originZ+(k*tau) );
		}
	}
}

