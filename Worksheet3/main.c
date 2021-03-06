#include "helper.h"
#include "visual.h"
#include "init.h"
#include <stdio.h>
#include"uvp.h"
#include"boundary_val.h"
#include"sor.c"
#include <sys/stat.h>
#include <dirent.h>

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args){

	double **U, **V, **P, **F, **G, **RS;
	int **Flag;
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
	Flag = imatrix(0, imax+1, 0, jmax+1);
	
	/* Assign initial values to u, v, p */
	init_uvp(UI, VI, PI, imax, jmax, U, V, P);
	
	/* Initialize lower part of the domain with UI = 0 for the flow_over_step problem */
	/* (this code might be moved to somewhere else later) */
	if(strcmp(problem, "flow_over_step") == 0) {
		init_matrix(U, 0, imax ,  0, jmax/2, 0);
	}

	/* Initialization of flag field */
	init_flag(problem, imax, jmax, Flag, dp);	

	clock_gettime(CLOCK_MONOTONIC, &currentTime);
	
	while(t <= t_end){
	
		/* Select δt */
		calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
		
		/* Set boundary values for u and v */
		boundaryvalues(imax, jmax, U, V, wl, wr, wt, wb, Flag);
	
		/* Set special boundary values */
		spec_boundary_val(problem, imax, jmax, U, V, P, Re, xlength, ylength, dp);

		/* Compute F(n) and G(n) */
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, Flag);
		
		/* Compute the right-hand side rs of the pressure equation */
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);
		
		/* Perform SOR iterations */
		it = 0;
		res = 1e6;
		while(it < itermax && res > eps){
			sor(omg, dx, dy, dp, imax, jmax, P, RS, &res, Flag);
			it++;
		}
		
		/* Compute u(n+1) and v(n+1) */
		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, Flag);
		
		t = t + dt;
		n++;
		
		/* Generate snapshot for current timestep */
		if((int) n % timestepsPerPlotting == 0) {
			write_vtkFile(output_dirname, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
		}
		
		/* Print out simulation time */
		printf("Time: %.4f", t);
		
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
