#include "helper.h"
#include "init.h"
#include <string.h>

int read_parameters( const char *szFileName,       /* name of the file */
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
		    int *wb)
{
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

/*
1. Supplement the function read parameters such that it reads initial values for the variable
problem and for the boundary value options wl,wr,wt and wb. (page 13, worksheet 3)
*/
   /*READ_STRING( szFileName, *problem);*/
   READ_INT (szFileName, *wl);
   READ_INT (szFileName, *wr);
   READ_INT (szFileName, *wt);
   READ_INT (szFileName, *wb);

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   return 1;
}


void init_uvp(
  double UI,
  double VI,
  double PI,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **P
){
/*
	init_matrix_U(U, 0, imax ,  0, jmax+1, UI);
	init_matrix_V(V, 0, imax+1, 0, jmax  , VI);
	init_matrix_P(P, 0, imax+1, 0, jmax+1, PI);
*/
}

void init_flag(
	char* problem, 
	int imax, 
	int jmax, 
	int **Flag,
	double dp
){
	char image_filename[60];
	int **pic;
	int i,j;
	
/*call read_pgm to this program,and check domain, and check the consistency of domain size with Flag matrix*/
	strcpy(image_filename, problem);
	strcat(image_filename, ".pgm");     
	pic=read_pgm(image_filename);	

/*	for(j=0;j<22;j++){
		for(i=0;i<22;i++){
		printf(" %d ",pic[i][j]);
		}
		printf("\n");
	}*/


	/*before this stage, pic setting: fluid cell is 1, obstacle is 0*/
	for(i=1;i<=imax;i++){
		for(j=1;j<=jmax;j++){
			Flag[i][j] = pic[i][j] * 16 + pic[i+1][j] * 8 + pic[i-1][j] * 4 + pic[i][j-1] * 2 + pic[i][j+1] * 1;
		}
	}
	
	if(strcmp(problem, "plane_shear_flow") == 0) {
		if(dp != 0){
			for(j = 1; j <= jmax; j++) {
				Flag[0][j] = Flag[0][j] | 64;
			}
		}
	}
	
/*	Here depending on "problem"
	for(j=1;j<=jmax;j++){
		Flag[0][j]=pic[0][j]*16+pic[0+1][j]*8+pic[0][j-1]*2+pic[0][j+1]*1;
	}
	for(j=1;j<=jmax;j++){
		Flag[imax+1][j]=pic[imax+1][j]*16+pic[imax+1-1][j]*4+pic[imax+1][j-1]*2+pic[imax+1][j+1]*1;
	}
	for(i=0;j<=imax+1;i++){
		Flag[i][0]=pic[i][0]*16+pic[i+1][0]*8+pic[i-1][0]*4+pic[i][0+1]*1;
	}
	for(i=0;j<=imax+1;i++){
		Flag[i][jmax+1]=pic[i][jmax+1]*16+pic[i+1][jmax+1]*8+pic[i-1][jmax+1]*4+pic[i][jmax+1-1]*2;
	}
*/


}

