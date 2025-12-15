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


void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int n = *la;
  double h = 1.0 / (n + 1);
  double factor = 1.0 / (h * h);

  for (int i = 0; i < n; i++){
    RHS[i] = 0.0;
  }
  if (n > 0){
    RHS[0]     = factor * (*BC0);
    RHS[n - 1] = factor * (*BC1);
  }
}  


void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la,double* BC0, double* BC1){
  int n = *la;
  double T0 = *BC0;
  double T1 = *BC1;

  for (int i = 0; i < n; i++){
    EX_SOL[i] = T0 + X[i] * (T1 - T0);
  }
}  


void set_grid_points_1D(double* x, int* la){
  int n = *la;
  double h = 1.0 / (n + 1);

  for (int i = 0; i < n; i++){
    x[i] = (i + 1) * h;
  }
}

double relative_forward_error(double* x, double* y, int* la){
  int n = *la;
  double num = 0.0;
  double den = 0.0;

  for (int i = 0; i < n; i++){
    double diff = x[i] - y[i];
    num += diff * diff;
    den += y[i] * y[i];
  }

  return sqrt(num / den);
}


int indexABCol(int row, int col, int *lab){
  int ldab = *lab;
  return row + col * ldab;
}


int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  int ncols = *n;
  int ldab = *lab;
  int diag = *kl + *ku;      /* Diagonal row index in band storage (0-based) */
  int sub  = diag + *kl;     /* Sub-diagonal row index */
  int super = diag - *ku;    /* Super-diagonal row index */

  *info = 0;
  if (ncols <= 0) {
    return *info;
  }

  for (int i = 0; i < ncols; i++) {
    ipiv[i] = i + 1; /* No pivoting for the tridiagonal case */
  }

  for (int j = 0; j < ncols - 1; j++) {
    double pivot = AB[j * ldab + diag];
    if (pivot == 0.0) {
      *info = j + 1; /* 1-based position of zero pivot */
      return *info;
    }
    double factor = AB[j * ldab + sub] / pivot;
    AB[j * ldab + sub] = factor;

    /* Update the next diagonal element */
    AB[(j + 1) * ldab + diag] -= factor * AB[(j + 1) * ldab + super];
  }

  if (AB[(ncols - 1) * ldab + diag] == 0.0) {
    *info = ncols;
  }
  return *info;
}