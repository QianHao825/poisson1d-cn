#include <math.h>   
#include <stdlib.h> 
#include <string.h>  
#include "lib_poisson1D.h"

/*********************/
/* EX10: CSR / CSC   */
/*********************/

static double l2_norm(int n, const double *v) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += v[i] * v[i];
  return sqrt(s);
}

int nnz_poisson1D(int n) {
  return (n <= 0) ? 0 : (3 * n - 2);
}

/* A = tridiag(-1, 2, -1) */
void set_CSR_poisson1D(int n, int *rowptr, int *colind, double *val) {
  int k = 0;
  rowptr[0] = 0;

  for (int i = 0; i < n; ++i) {
    /* left */
    if (i - 1 >= 0) {
      colind[k] = i - 1;
      val[k]    = -1.0;
      k++;
    }
    /* diag */
    colind[k] = i;
    val[k]    =  2.0;
    k++;

    /* right */
    if (i + 1 < n) {
      colind[k] = i + 1;
      val[k]    = -1.0;
      k++;
    }
    rowptr[i + 1] = k;
  }
}

void set_CSC_poisson1D(int n, int *colptr, int *rowind, double *val) {
  int k = 0;
  colptr[0] = 0;

  for (int j = 0; j < n; ++j) {
    /* up (row j-1, col j) */
    if (j - 1 >= 0) {
      rowind[k] = j - 1;
      val[k]    = -1.0;
      k++;
    }
    /* diag (row j, col j) */
    rowind[k] = j;
    val[k]    =  2.0;
    k++;

    /* down (row j+1, col j) */
    if (j + 1 < n) {
      rowind[k] = j + 1;
      val[k]    = -1.0;
      k++;
    }
    colptr[j + 1] = k;
  }
}

void dcsrmv(int n,
            const int *rowptr, const int *colind, const double *val,
            const double *x, double *y) {
  for (int i = 0; i < n; ++i) {
    double sum = 0.0;
    for (int k = rowptr[i]; k < rowptr[i + 1]; ++k) {
      sum += val[k] * x[colind[k]];
    }
    y[i] = sum;
  }
}

void dcscmv(int n,
            const int *colptr, const int *rowind, const double *val,
            const double *x, double *y) {
  /* y = 0 */
  for (int i = 0; i < n; ++i) y[i] = 0.0;

  for (int j = 0; j < n; ++j) {
    const double xj = x[j];
    for (int k = colptr[j]; k < colptr[j + 1]; ++k) {
      y[rowind[k]] += val[k] * xj;
    }
  }
}

/* Richardson: x <- x + alpha*(b - A x) */
void richardson_alpha_csr(int n, const int *rowptr, const int *colind, const double *val, const double *b, double *x, double alpha, double tol, int maxit, double *resvec, int *nbite) {
  double *Ax = (double*)malloc((size_t)n * sizeof(double));
  double *r  = (double*)malloc((size_t)n * sizeof(double));
  double nb  = l2_norm(n, b);
  if (nb == 0.0) nb = 1.0;

  int it = 0;
  for (; it < maxit; ++it) {
    dcsrmv(n, rowptr, colind, val, x, Ax);
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ax[i];

    double nr = l2_norm(n, r) / nb;
    resvec[it] = nr;
    if (nr < tol) break;

    for (int i = 0; i < n; ++i) x[i] += alpha * r[i];
  }

  *nbite = it;
  free(Ax);
  free(r);
}

void richardson_alpha_csc(int n, const int *colptr, const int *rowind, const double *val, const double *b, double *x, double alpha, double tol, int maxit, double *resvec, int *nbite) {
  double *Ax = (double*)malloc((size_t)n * sizeof(double));
  double *r  = (double*)malloc((size_t)n * sizeof(double));
  double nb  = l2_norm(n, b);
  if (nb == 0.0) nb = 1.0;

  int it = 0;
  for (; it < maxit; ++it) {
    dcscmv(n, colptr, rowind, val, x, Ax);
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ax[i];

    double nr = l2_norm(n, r) / nb;
    resvec[it] = nr;
    if (nr < tol) break;

    for (int i = 0; i < n; ++i) x[i] += alpha * r[i];
  }

  *nbite = it + 1;
  free(Ax);
  free(r);
}

/* Jacobi in CSR: x^{k+1} = D^{-1}(b - (L+U)x^k) */
void jacobi_csr(int n, const int *rowptr, const int *colind, const double *val, const double *b, double *x, double tol, int maxit, double *resvec, int *nbite) {
  double *xnew = (double*)malloc((size_t)n * sizeof(double));
  double *r    = (double*)malloc((size_t)n * sizeof(double));
  double nb    = l2_norm(n, b); if (nb == 0.0) nb = 1.0;

  int it = 0;
  for (; it < maxit; ++it) {
    for (int i = 0; i < n; ++i) {
      double diag = 0.0, sum = 0.0;
      for (int k = rowptr[i]; k < rowptr[i+1]; ++k) {
        int j = colind[k];
        if (j == i) diag = val[k];
        else sum += val[k] * x[j];
      }
      xnew[i] = (b[i] - sum) / diag;
    }

    /* residual r = b - A*xnew */
    dcsrmv(n, rowptr, colind, val, xnew, r);
    for (int i = 0; i < n; ++i) r[i] = b[i] - r[i];

    double nr = l2_norm(n, r) / nb;
    resvec[it] = nr;
    for (int i = 0; i < n; ++i) x[i] = xnew[i];
    if (nr < tol) break;
  }

  *nbite = it + 1;
  free(xnew);
  free(r);
}

/* Gauss-Seidel in CSR: solve (D+L)x^{k+1} = b - U x^k by forward sweep */
void gauss_seidel_csr(int n, const int *rowptr, const int *colind, const double *val, const double *b, double *x,double tol, int maxit, double *resvec, int *nbite) {
  double *Ax = (double*)malloc((size_t)n * sizeof(double));
  double *r  = (double*)malloc((size_t)n * sizeof(double));
  double nb  = l2_norm(n, b); if (nb == 0.0) nb = 1.0;

  int it = 0;
  for (; it < maxit; ++it) {
    for (int i = 0; i < n; ++i) {
      double diag = 0.0;
      double sumL = 0.0, sumU = 0.0;

      for (int k = rowptr[i]; k < rowptr[i+1]; ++k) {
        int j = colind[k];
        if (j == i) diag = val[k];
        else if (j < i) sumL += val[k] * x[j];  /* updated x */
        else            sumU += val[k] * x[j];  /* old x */
      }
      x[i] = (b[i] - sumL - sumU) / diag;
    }

    dcsrmv(n, rowptr, colind, val, x, Ax);
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ax[i];

    double nr = l2_norm(n, r) / nb;
    resvec[it] = nr;
    if (nr < tol) break;
  }

  *nbite = it + 1;
  free(Ax);
  free(r);
}
