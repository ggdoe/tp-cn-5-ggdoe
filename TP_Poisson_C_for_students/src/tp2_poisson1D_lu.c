/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int NRHS;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double *AB, *LB, *UB;
  double *A, *L, *U;
  double *y;
  double temp, relres;
  clock_t begin, end;
  const int nbr_rep = 1000;


  NRHS=1;
  // nbpoints=102;
  nbpoints= 10000;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  kv=0;// kv n'est pas requis pour dgbmv
  ku=1;// ku diagonale upper
  kl=1;// kl diagonale lower
  lab=kv+kl+ku+1; // ku + kl + diag principale


  printf("--------- LU BAND ---------\n");
  RHS=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  y = (double *) calloc(la, sizeof(double));

  set_grid_points_1D(X, &la); // X = discr de ]0,1[ avec 'la' points
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1); // RHS = (T0, 0, 0, ..., 0, T1)
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1); // 'y' est la solution exacte avec 'la' points
   


  //set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
// set_GB_operator_rowMajor_poisson1D(AB, &lab, &la, &kv);

  write_vec(RHS, &la, "RHS.dat");
  write_vec(y, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  // stockage BAND
  AB = (double *) malloc( lab  * la * sizeof(double)); // allocation matrice + set 0
  LB = (double *) malloc((1+kl)* la * sizeof(double)); // allocation matrice + set 0
  UB = (double *) malloc((1+ku)* la * sizeof(double)); // allocation matrice + set 0

  // myluB_rowMajor_poisson1D(LB, UB, &la);
  begin = clock();
  for(int i = 0; i < nbr_rep; i++)
    myluB_colMajor_poisson1D(LB, UB, &la);
  end = clock();

  // dgbmm n'existe pas ==> on calcule Uy puis L(Uy) et on compare à RHS
  cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, 0, ku, 1, UB, 1 + ku, EX_SOL, 1, 0, y, 1); // Ux = y
  cblas_dswap(la, EX_SOL, 1, y, 1); // on swap EX_SOL et 'y' pour plus de cohérence dans la notation dans la suite
  cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, 0, 1, LB, 1 + kl, EX_SOL, 1, 0, y, 1); // Lx =  y

  // cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, 0, ku, 1, UB, la, EX_SOL, 1, 0, y, 1); // Ux = y
  // cblas_dswap(la, EX_SOL, 1, y, 1); // on swap EX_SOL et 'y' pour plus de cohérence dans la notation dans la suite
  // cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, kl, 0, 1, LB, la, EX_SOL, 1, 0, y, 1); // Ly = y
  
  // for(int i = 0; i < la; i++)
  //   printf("%.0f ", y[i]);

//  printf("\n");
//   for(int i = 0; i < la; i++)
//     printf("%.4f ", RHS[i]);

  /* Relative residual */
  temp = cblas_ddot(la, RHS, 1, RHS,1);
  cblas_daxpy(la, -1.0, RHS, 1, y, 1);
  relres = cblas_ddot(la, y, 1, y,1);
  relres = sqrt(relres / temp);

  printf("\nThe relative residual error is relres = %e\n",relres);
  unsigned long micros = (end -  begin) * 10e6 / CLOCKS_PER_SEC / nbr_rep;
  printf("\n\nTime mesure = %ld µs\n",micros);

  printf("\n--------- LU DENSE ---------\n\n");
  // stockage DENSE
  A = (double *) calloc(la*la, sizeof(double)); // allocation matrice + set 0
  L = (double *) calloc(la*la, sizeof(double)); // allocation matrice + set 0
  U = (double *) calloc(la*la, sizeof(double)); // allocation matrice + set 0

  begin = clock();
  for(int i = 0; i < nbr_rep; i++)
    mylu_rowMajor_poisson1D(L, U, &la);
  end = clock();

  set_GE_rowMajor_operator_poisson1D(A, &la);

  // for(int i = 0; i < la * la; i++){
  //   if(i%la == 0)
  //     printf("\n");
  //   printf("%.0lf\t", U[i]);
  // }

  temp = LAPACKE_dlansy(LAPACK_ROW_MAJOR, 'F', 'U', la, A, la); // on profite de la symetrie de A
  // temp = LAPACKE_dlange(LAPACK_ROW_MAJOR, 'F', la, la, A, la);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, la, la, la, 1, L, la, U, la,-1, A, la);
  relres = LAPACKE_dlange(LAPACK_ROW_MAJOR, 'F', la, la, A, la);
  relres = relres/temp;

  printf("The relative residual error is relres = %e",relres);
  micros = (end -  begin) * 10e6 / CLOCKS_PER_SEC / nbr_rep;
  printf("\n\nTime mesure = %ld µs\n",micros);
// LAPACKE_dlansy // SY

// cblas_dswap()
  // cblas_dgemm(); // LU-A
  // cblas_dgemm(); // ||LU-A||
  // cblas_dgemm(); // ||A||
  
  

  free(RHS);
  free(EX_SOL);
  free(X);
  free(y);
  free(AB); free(UB); free(LB);
  free(A);  free(U);  free(L);

  printf("\n\n--------- End -----------\n");
}
