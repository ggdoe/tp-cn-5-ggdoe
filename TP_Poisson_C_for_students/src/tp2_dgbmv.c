/******************************************/
/* tp2_dgbmv.c                 */
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
  double *AB;
  double *y;
  double temp, relres;
  clock_t begin, end;
  const int nbr_rep = 10000;

  NRHS=1;
  // nbpoints=102;
  nbpoints=10000;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  kv=0;// kv n'est pas requis pour dgbmv
  ku=1;// ku diagonale upper
  kl=1;// kl diagonale lower
  lab=kv+kl+ku+1; // ku + kl + diag principale


  printf("--------- DGBMV ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);
  y = (double *) calloc(la, sizeof(double));


  kv = 0; // dgbmv n'effectue pas de factorisation LU, il n'y a pas de 0 de padding
  set_grid_points_1D(X, &la); // X = discr de ]0,1[ avec 'la' points
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1); // RHS = (T0, 0, 0, ..., 0, T1)
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1); // EX_SOL solution exacte avec 'la' points
   
   
  AB = (double *) malloc(sizeof(double)*lab*la); // allocation matrice

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
// set_GB_operator_rowMajor_poisson1D(AB, &lab, &la, &kv);

  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  begin = clock(); 
  for(int i = 0; i < nbr_rep; i++)
    cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1, AB, lab, EX_SOL, 1, 0, y, 1);// y est nul
  end = clock(); 

  // rowMajor ne fonctionne tjr pas
  // cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, kl, ku, 1, AB, la, EX_SOL, 1, 0, y, 1);
  
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
////////////////////// end ex.4


  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(y);

  printf("\n\n--------- End -----------\n");
}
