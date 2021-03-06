#include "helper.h"
#include "init.h"
#include <string.h>
#include <mpi.h>
#include "boundary_val.h"

int read_parameters( const char *szFileName,       /* name of the file */
					char* pgm,      /* name of pgm file */
                    double *Re,                /* reynolds number   */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dt,                /* time step */
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations  */
		                               /* for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
		    double *dt_value,           /* time for output */
		    char *problem,
		    double *dp,
	            int *wl,
		    int *wr,
		    int *wt,
		    int *wb,
		    int *timestepsPerPlotting,
		    int *iproc,
		    int *jproc)
{

   READ_STRING( szFileName, pgm );

   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );

   READ_DOUBLE( szFileName, *Re    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *dt    );

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *alpha );

   READ_INT   ( szFileName, *itermax );
   READ_DOUBLE( szFileName, *dt_value );

   READ_DOUBLE( szFileName, *dp );
   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *PI );

   READ_INT (szFileName, *wl);
   READ_INT (szFileName, *wr);
   READ_INT (szFileName, *wt);
   READ_INT (szFileName, *wb);
   
   READ_INT (szFileName, *timestepsPerPlotting);
   
   READ_INT (szFileName, *iproc);
   READ_INT (szFileName, *jproc);

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   return 1;
}


void init_uvp(
  double UI,
  double VI,
  double PI,
  int il, int ir,
  int jb, int jt,
  double **U,
  double **V,
  double **P
){
	init_matrix(U, il-2, ir+1, jb-1, jt+1, UI);
	init_matrix(V, il-1, ir+1, jb-2, jt+1, VI);
	init_matrix(P, il-1, ir+1, jb-1, jt+1, PI);
}

void init_flag(
	char* pgm, 
	int imax,
	int jmax,
	int il, int ir,
	int jb, int jt,
	int **Flag,
	double dp
){
	int **pic;
	int i,j;
	
	if(strcmp(pgm, "-") == 0) {
	
		init_imatrix(Flag, il-1, ir+1, jb-1, jt+1, C_F);
		
	} else {
	
		pic = read_pgm(pgm);	

		/* Picture values: 1 for fluid, 0 for obstacle */
		for(i = max(il-1, 1); i <= min(ir+1, imax); i++){
			for(j =  max(jb-1, 1); j <= min(jt+1, jmax); j++){
				Flag[i][j] = max(pic[i+1][j] * B_O + pic[i-1][j] * B_W + pic[i][j-1] * B_S + pic[i][j+1] * B_N, pic[i][j] * C_F);
			}
		}
		
	}
	
	/* Pressure boundary conditions are currently implemented in such a way 
		that we don't need extra flags for them (might be implemented later) */
	/* if(dp != 0){
		for(j = 1; j <= jmax; j++) {
			Flag[0][j] = Flag[0][j] | 64;
		}
	} */

}

