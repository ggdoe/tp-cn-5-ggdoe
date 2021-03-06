/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_rowMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii;

  const int v = (*kv) * (*la); // nbr lignes laissée vide
  const int n = *la; // les premiers elements ne sont pas set

  // on n'a pas besoin d'initialiser les premiers éléments *
  for(ii = 1; ii < n; ii++){
    AB[v + ii] = -1;
  }
  for(ii = 0; ii < n; ii++){
    AB[v + n + ii] = 2;
  }
  for(ii = 0; ii < n-1; ii++){
    AB[v + 2*n + ii] = -1;
  }
  // on n'a pas besoin d'initialiser les derniers éléments *
}

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}
  
  AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii, jj;
  file = fopen(filename, "w");

  if(file != NULL){
      // TODO ...
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void eig_poisson1D(double* eigval, int *la){
  int ii;
  double scal;
  for (ii=0; ii< *la; ii++){
    scal=(1.0*ii+1.0)*M_PI_2*(1.0/(*la+1));
    eigval[ii]=sin(scal);
    eigval[ii]=4*eigval[ii]*eigval[ii];
  } 
}

double eigmax_poisson1D(int *la){
  double eigmax;
  eigmax=sin(*la *M_PI_2*(1.0/(*la+1)));
  eigmax=4*eigmax*eigmax;
  return eigmax;
}

double eigmin_poisson1D(int *la){
  double eigmin;
  eigmin=sin(M_PI_2*(1.0/(*la+1)));
  eigmin=4*eigmin*eigmin;
  return eigmin;
}

double richardson_alpha_opt(int *la){
  //TODO
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit){
  //TODO
}

void myluB_rowMajor_poisson1D(double *LB, double *UB, int *la){
  LB[0] = 1;
  for(int k = 1; k < *la; k++){
    UB[k] = -1;
    LB[k] = 1;
  }
  UB[*la] = 2;
  for(int k = 0; k < *la - 1; k++){
    const double inv_un = 1/UB[*la + k];
    LB[*la + k] = -inv_un;
    UB[*la + k + 1] = 2 - inv_un;
  }
}

void myluB_colMajor_poisson1D(double *LB, double *UB, int *la){
  UB[1] = 2;   // on n'initialise pas les *
  for(int k = 1; k < *la; k++){
    UB[2*k] = -1;
    LB[2*(k-1)] = 1;
    const double inv_un = 1/UB[2*k-1];
    UB[2*k+1] = 2 - inv_un;
    LB[2*k-1] = -inv_un; // on décale les indices pour profiter de inv_un
  }
  LB[2 * (*la-1)] = 1;
}

void mylu_rowMajor_poisson1D(double *L, double *U, int *la){
  // L et U doivent etre initialiser à 0 avec calloc.
  for(int i = 0; i < *la; i++){
    L[i * (*la) + i] = 1;
  }
  for(int i = 0; i < *la-1; i++){
    U[i * (*la) + i + 1] = -1;
  }

  U[0] = 2;
  for(int i = 0; i < *la - 1; i++){
    const double inv_un = 1/U[i * (*la) + i];
    U[(i+1) * (*la) + i + 1] = 2 - inv_un;
    L[(i+1) * (*la) + i] = -inv_un;
  }
}

void set_GE_rowMajor_operator_poisson1D(double* A, int *la){
  // A doit etre initialiser à 0 avec calloc.
  A[0] = 2;
  for(int i = 1; i < *la; i++){
    A[(i-1) * (*la) + i] = -1;
    A[i * (*la) + i] = 2;
    A[i * (*la) + i - 1] = -1;
  }
}
