/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int n    = *la;         
  int ldab = *lab;         

  double h    = 1.0 / (n + 1);
  double diag =  2.0 / (h * h);
  double off  = -1.0 / (h * h);

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < ldab; i++) {
      AB[indexABCol(i, j, lab)] = 0.0;
    }
  }

  int row_extra  = 0; 
  (void)row_extra;   
  int row_super  = 1;  
  int row_diag   = 2;  
  int row_sub    = 3;  

  for (int j = 0; j < n; j++) {
    if (j > 0) {
      AB[indexABCol(row_super, j, lab)] = off;
    }

    AB[indexABCol(row_diag, j, lab)] = diag;

    if (j < n - 1) {
      AB[indexABCol(row_sub, j, lab)] = off;
    }
  }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int n    = *la;
  int ldab = *lab;


  for (int j = 0; j < n; j++) {
    for (int i = 0; i < ldab; i++) {
      AB[indexABCol(i, j, lab)] = 0.0;
    }
  }

  int row_diag = 2;

  for (int j = 0; j < n; j++) {
    AB[indexABCol(row_diag, j, lab)] = 1.0;
  }
}


void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the identity matrix
  // Only the main diagonal should have 1, all other entries are 0
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // TODO: Compute RHS vector
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D
}  

void set_grid_points_1D(double* x, int* la){
  // TODO: Generate uniformly spaced grid points in [0,1]
}

double relative_forward_error(double* x, double* y, int* la){
  // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)
  return 0.0;
}

int indexABCol(int i, int j, int *lab){
  // TODO: Return the correct index formula for column-major band storage
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices
  return *info;
}