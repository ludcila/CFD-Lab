#ifndef __SOR_H_
#define __SOR_H_

/**
 * One GS iteration for the pressure Poisson equation. Besides, the routine must 
 * also set the boundary values for P according to the specification. The 
 * residual for the termination criteria has to be stored in res.
 * 
 * An \omega = 1 GS - implementation is given within sor.c.
 */
void sor(
  double omg,
  double dx,
  double dy,
  double dp,
  int il, int ir,
  int jb, int jt,
  int    imax,
  int    jmax,
  int rank_l, int rank_r, 
  int rank_b, int rank_t,
  double **P,
  double **RS,
  double *res,
  int **Flag,
  double *bufSend,
  double *bufRecv
);


#endif
