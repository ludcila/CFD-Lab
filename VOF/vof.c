/* Set the initial fluid fraction based on the pgm file */
void init_fluidFraction(int **pgm, double **fluidFraction, int imax, int jmax) {

}

/* Bookeeping adjustments as described by Hirt and Nichols */
void adjust_fluidFraction(double **fluidFraction, int **flagField, double epsilon, int imax, int jmax) {

}

/* Determine orientation of the free surface */
void calculate_freeSurfaceOrientation(double **fluidFraction, double **dFdx, double **dFdy, int imax, int jmax) {

}

/* Timestepping for the fluid fraction field */
void calculate_fluidFraction(
	double **fluidFraction, 
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

}
