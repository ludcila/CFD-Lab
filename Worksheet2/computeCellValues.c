#include "computeCellValues.h"
#include "LBDefinitions.h"

void computeDensity(const double *const currentCell, double *density){

	int i;
	
	for(i = 0; i < Q; i++) {
		*density += currentCell[i];
	}

}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){

	int i;
	
	for(i = 0; i < Q; i++) {
		velocity[0] += currentCell[i] * LATTICEVELOCITIES[i][0];
		velocity[1] += currentCell[i] * LATTICEVELOCITIES[i][1];
		velocity[2] += currentCell[i] * LATTICEVELOCITIES[i][2];
	}

}

void computeFeq(const double * const density, const double * const velocity, double *feq){

	int i;
	double u_dot_u = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
	
	for(i = 0; i < Q; i++) {
		double c_dot_u = LATTICEVELOCITIES[i][0] * velocity[0] + LATTICEVELOCITIES[i][1] * velocity[1] + LATTICEVELOCITIES[i][2] * velocity[2];
		double T1 = c_dot_u / pow(C_S, 2.0);
		double T2 = pow(c_dot_u, 2.0) / 2.0 / pow(C_S, 4.0);
		double T3 = u_dot_u / 2.0 / pow(C_S, 2.0);
		feq[i] = LATTICEWEIGHTS[i] * *density * (1 + T1 + T2 - T3);
	}

}

