#ifndef __VOF__

#define __VOF__

/* Flags for free surface cells */
enum FS_xy {
	FS_N 	= 32+16+1,		/* Empty cell at north */
	FS_S	= 32+16+2,		/* Empty cell at south */
	FS_O	= 32+16+4,		/* Empty cell at east */
	FS_W	= 32+16+8		/* Empty cell at west */
};
/* Flags to indicate cell type */
enum C_type {
	C_F 	= 16,		/* Interior fluid cell */
	C_FS	= 48,		/* Free surface cell */
	C_E		= 0			/* Empty cell */
};

#endif

/* Set the initial fluid fraction based on the pgm file */
void init_fluidFraction(int **pgm, double **fluidFraction, int imax, int jmax);

/* Bookeeping adjustments as described by Hirt and Nichols */
void adjust_fluidFraction(double **fluidFraction, int **flagField, double epsilon, int imax, int jmax);

/* Determine orientation of the free surface */
void calculate_freeSurfaceOrientation(double **fluidFraction, double **dFdx, double **dFdy, int imax, int jmax);

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
);
