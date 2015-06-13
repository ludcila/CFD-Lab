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
