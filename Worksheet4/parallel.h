#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

#define		LEFT_TO_RIGHT	1
#define		RIGHT_TO_LEFT	2
#define		UP_TO_DOWN		3
#define		DOWN_TO_UP		4

#define		VAR_P			1
#define		VAR_U			2
#define		VAR_V			3

int get_location_rank(int omg_i, int omg_j, int iproc, int jproc);

void init_parallel(
	/* Global domain size */
	int iproc, 
	int jproc, 
	int imax, 
	int jmax,
	/* Rank and local index ranges */
	int *myrank,
	int *il,
	int *ir,
	int *jb,
	int *jt,
	/* Rank of neighbours */
	int *rank_l,
	int *rank_r,
	int *rank_b,
	int *rank_t,
	/* Global block index (?) */
	int *omg_i,
	int *omg_j,
	int num_proc
);


void Program_Message(char *txt);
/* produces a stderr text output  */



void Programm_Sync(char *txt);
/* produces a stderr textoutput and synchronize all processes */



void Programm_Stop(char *txt);
/* all processes will produce a text output, be synchronized and finished */

void copy_to_buffer(
	double **matrix, 
	double *buffer, 
	int il, int ir, 
	int jb, int jt
);

void copy_from_buffer(
	double **matrix,
	double *buffer,
	int il, int ir,
	int jb, int jt
);

void exchange(
	double **matrix,
	int il, int ir,
	int jb, int jt,
	int rank_l, int rank_r,
	int rank_b, int rank_t,
	int direction,
	int variable,
	double *bufSend, double *bufRecv
);

void pressure_com(
	double **P, 
	int il, int ir, 
	int jb, int jt, 
	int rank_l, int rank_r, 
	int rank_b, int rank_t,
	double *bufSend, double *bufRecv
);

void uv_com(
	double **U, 
	double **V, 
	int il, int ir, 
	int jb, int jt, 
	int rank_l, int rank_r,
	int rank_b, int rank_t,
	double *bufSend, double *bufRecv
);
	
