/**********************************************/
/* lib_poisson1D.h                            */
/* Header for Numerical library developed to  */ 
/* solve 1D Poisson problem (Heat equation)   */
/**********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "blaslapack_headers.h"

void set_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int *la);
void set_GB_operator_colMajor_poisson1D(double* AB, int* lab, int *la, int *kv);
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int* lab, int *la, int *kv);
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1);
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1);
void set_grid_points_1D(double* x, int* la);
void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int *la, char* filename);
void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename);
void write_vec(double* vec, int* la, char* filename);
void write_xy(double* vec, double* x, int* la, char* filename);
void eig_poisson1D(double* eigval, int *la);
double eigmax_poisson1D(int *la);
double eigmin_poisson1D(int *la);
double richardson_alpha_opt(int *la);
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit);
