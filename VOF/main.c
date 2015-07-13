#include "helper.h"
#include "visual.h"
#include "init.h"
#include <stdio.h>
#include"uvp.h"
#include"boundary_val.h"
#include"sor.c"
#include "vof.h"
#include <dirent.h>
#include <sys/stat.h>

int main(int argn, char** args){

	double **U, **V, **P, **F, **G, **RS;
	const char *szFileName = "dam_break.dat";
	double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value;
	double res = 0, t = 0, n = 0;
	int imax, jmax, itermax, it;

	
	/* Additional data structures for VOF */
	double **fluidFraction;
	double **fluidFraction_alt;
	int **flagField;
	double **dFdx, **dFdy;
	int **pic;
	char output_dirname[60];
	double epsilon = 1e-3;
	
	/* Read the program configuration file using read_parameters() */
	read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value);        

	/* Set up the matrices (arrays) needed using the matrix() command */
	U = matrix(0, imax  , 0, jmax+1);
	V = matrix(0, imax+1, 0, jmax  );
	P = matrix(0, imax+1, 0, jmax+1);
	F = matrix(0, imax  , 0, jmax+1);
	G = matrix(0, imax+1, 0, jmax  );
	RS= matrix(0, imax+1, 0, jmax+1);
	flagField = imatrix(0, imax+1, 0, jmax+1);
	fluidFraction = matrix(0, imax+1, 0, jmax+1);
	fluidFraction_alt = matrix(0, imax+1, 0, jmax+1);
	dFdx = matrix(0, imax+1, 0, jmax+1);
	dFdy = matrix(0, imax+1, 0, jmax+1);


	/*create a directory*/
	strcpy(output_dirname, "dam_break/dam_break");
	mkdir("dam_break", 0777);
	
	/* Read pgm file with domain and initial setting */
	pic = read_pgm("dam_break.pgm");
	init_fluidFraction(pic, fluidFraction,fluidFraction_alt, imax, jmax);
	
	/* Assign initial values to u, v, p */
	init_uvp(UI, VI, PI, imax, jmax, U, V, P);

		/* Set valid values for the fluid fraction */
		adjust_fluidFraction(fluidFraction, flagField, epsilon, imax, jmax);

	while(t <= t_end){   
	
		/*Select δt*/
		calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
		
		/* Set boundary values for u and v */
/*		boundaryvalues(imax, jmax, U, V, flagField, dx, dy);
*/		
		/* Compute F(n) and G(n) */
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, flagField);
		
		/* Compute the right-hand side rs of the pressure equation */
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS, flagField);



		
				
		
		/* Determine the orientation of the free surfaces */

		calculate_freeSurfaceOrientation(fluidFraction, flagField, dFdx, dFdy, dx, dy, imax, jmax);

		
		/* Perform SOR iterations */
		it=0;
		res = 1e6;
		while(it < itermax && res > eps){
			sor(omg, dx, dy, imax, jmax, P, RS, fluidFraction, flagField, dFdx, dFdy, &res);
			it++;
		}

		
		/* Compute u(n+1) and v(n+1) */
		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, flagField);

		boundaryvalues(imax, jmax, U, V, flagField, dx, dy);
		
		/* Compute fluidFraction(n+1) */
		calculate_fluidFraction(fluidFraction,fluidFraction_alt, U, V, dFdx, dFdy, imax, jmax, dx, dy, dt);


		/* Set valid values for the fluid fraction */
		adjust_fluidFraction(fluidFraction, flagField, epsilon, imax, jmax);



		if((int)n % 20 == 0) {
			write_vtkFile(output_dirname, n, xlength, ylength, imax, jmax, dx, dy, U, V, P, fluidFraction, flagField);
		}
		
		/* Print out simulation time and whether SOR converged */
		printf("Time: %.4f", t);
		if(res > eps) {printf("\t*** Did not converge (res=%f, eps=%f)", res, eps);return -1;}
		printf("\n");
		
		t = t + dt;
		n++;
		if(t>0.9)
			printf("time > 0.9");
	
	}
	

	free_matrix(U , 0, imax  , 0, jmax+1);
	free_matrix(V , 0, imax+1, 0, jmax  );
	free_matrix(P , 0, imax+1, 0, jmax+1);
	free_matrix(F , 0, imax  , 0, jmax+1);
	free_matrix(G , 0, imax+1, 0, jmax  );
	free_matrix(RS, 0, imax+1, 0, jmax+1);
	
	return -1;
	
}
