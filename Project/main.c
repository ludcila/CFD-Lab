#include "helper.h"
#include "visual.h"
#include "init.h"
#include <stdio.h>
#include"uvp.h"
#include"boundary_val.h"
#include"sor.c"
#include <sys/stat.h>
#include <dirent.h>
#include "vof.h"

int main(int argn, char** args){

	double **U, **V, **P, **F, **G, **RS;
	int **flagField;
	char problem[60];
	char parameters_filename[60];
	char output_dirname[60];
	double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value, dp;
	double res = 0, t = 0, n = 0;
	int imax, jmax, itermax, it;
	int wl, wr, wt, wb;
	int timestepsPerPlotting;
	char old_output_filename[128];
	struct dirent *old_outputfile;
	DIR *output_dir;
	struct timespec previousTime, currentTime;
	double totalTime = 0;
	
	/* Additional data structures for VOF */
	double **fluidFraction;
	double **fluidFraction_alt;
	double **dFdx, **dFdy;
	double epsilon = 1e-10;

	/* Read name of the problem from the command line arguments */
	if(argn > 1) {
		strcpy(problem, args[1]);
	} else {
		printf("*** Please provide the name of the problem\n*** e.g. Run ./sim problem_name if there is a problem_name.dat file.\n");
		return 1;
	}

	/* Generate input filename based on problem name */
	strcpy(parameters_filename, problem);
	strcat(parameters_filename, ".dat");

	/* Read the program configuration file using read_parameters() */
	read_parameters(parameters_filename, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, problem, &dp, &wl, &wr, &wt, &wb, &timestepsPerPlotting);

	/* Create folder with the name of the problem */
	strcpy(output_dirname, problem);
	strcat(output_dirname, "/");
	strcat(output_dirname, problem);
	mkdir(problem, 0777);
	output_dir = opendir(problem);
	
	/* Delete existing files in output folder*/
	while((old_outputfile = readdir(output_dir))) {
		sprintf(old_output_filename, "%s/%s", problem, old_outputfile->d_name);
		remove(old_output_filename);
	}

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
	
	/* Assign initial values to u, v, p */
	init_uvp(UI, VI, PI, imax, jmax, U, V, P);
	
	/* Initialize lower part of the domain with UI = 0 for the flow_over_step problem */
	/* (this code might be moved to somewhere else later) */
	if(strcmp(problem, "flow_over_step") == 0) {
		init_matrix(U, 0, imax ,  0, jmax/2, 0);
	}

	/* Initialization of flag field */
	init_flag(problem, imax, jmax, flagField, dp);	

	/* Initialization of fluid fraction */
	init_fluidFraction(flagField, fluidFraction, fluidFraction_alt, imax, jmax);
	
	/* Set valid values for the fluid fraction */
	adjust_fluidFraction(fluidFraction, flagField, epsilon, imax, jmax);

	clock_gettime(CLOCK_MONOTONIC, &currentTime);
	
	while(t <= t_end){
	
		/* Generate snapshot for current timestep */
		if((int) n % timestepsPerPlotting == 0) {
			write_vtkFile(output_dirname, n, xlength, ylength, imax, jmax, dx, dy, U, V, P, fluidFraction, flagField);
		}
		
		/* Select δt */
		calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);

		/* Compute F(n) and G(n) */
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, flagField);
		
		/* Compute the right-hand side rs of the pressure equation */
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);
		
		/* Determine the orientation of the free surfaces */
		calculate_freeSurfaceOrientation(fluidFraction, flagField, dFdx, dFdy, dx, dy, imax, jmax);
		
		/* Perform SOR iterations */
		it = 0;
		res = 1e6;
		while(it < itermax && res > eps){
			sor(omg, dx, dy, dp, imax, jmax, P, RS, fluidFraction, flagField, dFdx, dFdy, &res);
			it++;
		}
		
		/* Compute u(n+1) and v(n+1) */
		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, flagField);
		
		/* Set boundary values for u and v */
		boundaryvalues(imax, jmax, U, V, wl, wr, wt, wb, flagField, dx, dy);
		
		/* TODO: Special boundary values? Do we use them with free surfaces? */
		
		/* Compute fluidFraction(n+1) */
		calculate_fluidFraction(fluidFraction,fluidFraction_alt, flagField, U, V, dFdx, dFdy, imax, jmax, dx, dy, dt);
		
		/* Set boundary values for u and v */
		boundaryvalues(imax, jmax, U, V, wl, wr, wt, wb, flagField, dx, dy);

		/* Set valid values for the fluid fraction */
		adjust_fluidFraction(fluidFraction, flagField, epsilon, imax, jmax);
		
		/* Set boundary values for u and v */
		boundaryvalues(imax, jmax, U, V, wl, wr, wt, wb, flagField, dx, dy);
		
		t = t + dt;
		n++;
		
		/* Print out simulation time */
		printf("(%d) Time: %.4f", (int)n, t);
		
		/* Print out runtime */
		previousTime = currentTime;
		clock_gettime(CLOCK_MONOTONIC, &currentTime);
		totalTime += (double)currentTime.tv_sec + 1e-9 * currentTime.tv_nsec - (double)previousTime.tv_sec - 1e-9 * previousTime.tv_nsec;
		printf("\tRuntime: %.3f s (avg runtime/step: %.3f s)", totalTime, totalTime/n);
		
		/* Print out whether the SOR converged */
		if(res > eps) printf("\t*** Did not converge (res=%f, eps=%f)", res, eps);
		printf("\n");
	
	}
	
	/* Close the output folder */
	closedir(output_dir);
	
	/* Tell user where to find the output */
	printf("Please find the output in the folder \"%s\".\n", problem);
	
	/* Free allocated memory */
	free_matrix(U , 0, imax  , 0, jmax+1);
	free_matrix(V , 0, imax+1, 0, jmax  );
	free_matrix(P , 0, imax+1, 0, jmax+1);
	free_matrix(F , 0, imax  , 0, jmax+1);
	free_matrix(G , 0, imax+1, 0, jmax  );
	free_matrix(RS, 0, imax+1, 0, jmax+1);
	
	return -1;
	
}
