#include "parallel.h"

int get_location_rank(int omg_i, int omg_j, int iproc, int jproc) {
	int proc_num = MPI_PROC_NULL;
	if(omg_i >= 0 && omg_i < iproc && omg_j >= 0 && omg_j < jproc) {
		proc_num = omg_i + omg_j * iproc;
	}
	return proc_num;
}

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
) {

	int subdomain_width = imax / iproc;
	int subdomain_height = jmax / jproc;
	int extra_width = imax % iproc;
	int extra_height = jmax % jproc;
	
	/* Get rank of this process */
	MPI_Comm_rank(MPI_COMM_WORLD, myrank);
	
	/* Get location of this subdomain */
	*omg_i = *myrank % iproc;
	*omg_j = *myrank / iproc;
	
	/* Get rank of neighbours */
	*rank_l = get_location_rank(*omg_i - 1, *omg_j, iproc, jproc);
	*rank_r = get_location_rank(*omg_i + 1, *omg_j, iproc, jproc);
	*rank_b = get_location_rank(*omg_i, *omg_j - 1, iproc, jproc);
	*rank_t = get_location_rank(*omg_i, *omg_j + 1, iproc, jproc);
	
	/* Get subdomain indices for this process */
	*il = *omg_i * subdomain_width + 1;
	*ir = *il + subdomain_width - 1;
	*jb = *omg_j * subdomain_height + 1;
	*jt = *jb + subdomain_height - 1;
	if(*rank_r == MPI_PROC_NULL) {
		*ir = *ir + extra_width;
	}
	if(*rank_t == MPI_PROC_NULL) {
		*jt = *jt + extra_height;
	}
	
	printf("myrank=%d location=(%d, %d) L=%d R=%d B=%d T=%d subdomain=[%d, %d]x[%d, %d]\n ", *myrank, *omg_i, *omg_j, *rank_l, *rank_r, *rank_b, *rank_t, *il, *ir, *jb, *jt);

}

void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}

void copy_to_strip(double **matrix, double *strip, int il, int ir, int jb, int jt) {
	int i, j;
	int k = 0;
	for(i = il; i <= ir; i++) {
		for(j = jb; j <= jt; j++) {
			strip[k++] = matrix[i][j];
		}
	}
}

void copy_from_strip(double **matrix, double *strip, int il, int ir, int jb, int jt) {
	int i, j;
	int k = 0;
	for(i = il; i <= ir; i++) {
		for(j = jb; j <= jt; j++) {
			matrix[i][j] = strip[k++];
		}
	}
}

void exchange(double **matrix, int il, int ir, int jb, int jt, int rank_l, int rank_r, int rank_b, int rank_t, int direction) {

	int num_elem;
	double *strip;
	int send_to, receive_from;
	int send_i_low, send_i_high, send_j_low, send_j_high;
	int receive_i_low, receive_i_high, receive_j_low, receive_j_high;
	
	switch(direction) {
		case LEFT_TO_RIGHT:
			num_elem = jt - jb + 1;
			send_to = rank_r;
			receive_from = rank_l;
			send_i_low = ir;
			send_i_high = ir;
			receive_i_low = il-1;
			receive_i_high = il-1;
			send_j_low = jb;
			send_j_high = jt;
			receive_j_low = jb;
			receive_j_high = jt;
			break;
		case RIGHT_TO_LEFT:
			num_elem = jt - jb + 1;
			send_to = rank_l;
			receive_from = rank_r;
			send_i_low = il;
			send_i_high = il;
			receive_i_low = ir+1;
			receive_i_high = ir+1;
			send_j_low = jb;
			send_j_high = jt;
			receive_j_low = jb;
			receive_j_high = jt;
			break;
		case DOWN_TO_UP:
			num_elem = ir - il + 1;
			send_to = rank_t;
			receive_from = rank_b;
			send_i_low = il;
			send_i_high = ir;
			receive_i_low = il;
			receive_i_high = ir;
			send_j_low = jt;
			send_j_high = jt;
			receive_j_low = jb-1;
			receive_j_high = jb-1;
			break;
		case UP_TO_DOWN:
			num_elem = ir - il + 1;
			send_to = rank_b;
			receive_from = rank_t;
			send_i_low = il;
			send_i_high = ir;
			receive_i_low = il;
			receive_i_high = ir;
			send_j_low = jb;
			send_j_high = jb;
			receive_j_low = jt+1;
			receive_j_high = jt+1;
			break;
	}
	strip = malloc(num_elem * sizeof(double));
	copy_to_strip(matrix, strip, send_i_low, send_i_high, send_j_low, send_j_high);
	MPI_Send(strip, num_elem, MPI_DOUBLE, send_to, 0, MPI_COMM_WORLD);
	MPI_Recv(strip, num_elem, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	copy_from_strip(matrix, strip, receive_i_low, receive_i_high, receive_j_low, receive_j_high);	
	free(strip);
	
}
