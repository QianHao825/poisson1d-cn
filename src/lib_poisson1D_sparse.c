#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "lib_poisson1D.h"

/* EX10: CSR / CSC   */

static double l2_norm(int n, const double *v) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += v[i] * v[i];
  return sqrt(s);
}

int nnz_poisson1D(int n) {
  /* tridiagonal: 3n-2 */
  if (n <= 1) return n;
  return 3 * n - 2;
}

/* Build Poisson 1D in CSR with coeffs: diag=2, off=-1 (NO 1/h^2 here) */
void set_CSR_poisson1D(int n, int *rowptr, int *colind, double *val) {
  int k = 0;
  rowptr[0] = 0;

  for (int i = 0; i < n; ++i) {
    if (i - 1 >= 0) {
      colind[k] = i - 1;
      val[k]    = -1.0;
      k++;
    }

    colind[k] = i;
    val[k]    = 2.0;
    k++;

    if (i + 1 < n) {
      colind[k] = i + 1;
      val[k]    = -1.0;
      k++;
    }
    rowptr[i + 1] = k;
  }
}

/* Build Poisson 1D in CSC with coeffs: diag=2, off=-1 (NO 1/h^2 here) */
void set_CSC_poisson1D(int n, int *colptr, int *rowind, double *val) {
  int k = 0;
  colptr[0] = 0;

  for (int j = 0; j < n; ++j) {
    if (j - 1 >= 0) {
      rowind[k] = j - 1;
      val[k]    = -1.0;
      k++;
    }

    rowind[k] = j;
    val[k]    = 2.0;
    k++;

    if (j + 1 < n) {
      rowind[k] = j + 1;
      val[k]    = -1.0;
      k++;
    }
    colptr[j + 1] = k;
  }
}

/* y = A*x (CSR) */
void dcsrmv(int n, const int *rowptr, const int *colind,
            const double *val, const double *x, double *y) {
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int k = rowptr[i]; k < rowptr[i+1]; ++k) {
      s += val[k] * x[colind[k]];
    }
    y[i] = s;
  }
}

/* y = A*x (CSC) */
void dcscmv(int n, const int *colptr, const int *rowind,
            const double *val, const double *x, double *y) {
  for (int i = 0; i < n; ++i) y[i] = 0.0;

  for (int j = 0; j < n; ++j) {
    double xj = x[j];
    for (int k = colptr[j]; k < colptr[j+1]; ++k) {
      y[rowind[k]] += val[k] * xj;
    }
  }
}

/* Richardson (CSR): x^{k+1} = x^k + alpha * (b - A x^k) */
void richardson_alpha_csr(int n, const int *rowptr, const int *colind,
                          const double *val, const double *b, double *x,
                          double alpha, double tol, int maxit,
                          double *resvec, int *nbite) {
  double *Ax = (double*)malloc((size_t)n * sizeof(double));
  double *r  = (double*)malloc((size_t)n * sizeof(double));

  double nb = l2_norm(n, b);
  if (nb == 0.0) nb = 1.0;

  int it = 0;
  for (; it < maxit; ++it) {
    dcsrmv(n, rowptr, colind, val, x, Ax);
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ax[i];

    double nr = l2_norm(n, r) / nb;
    resvec[it] = nr;
    if (nr < tol) { it++; break; } 

    for (int i = 0; i < n; ++i) x[i] += alpha * r[i];
  }

  *nbite = it;
  free(Ax);
  free(r);
}


void richardson_alpha_csc(int n, const int *colptr, const int *rowind,
                          const double *val, const double *b, double *x,
                          double alpha, double tol, int maxit,
                          double *resvec, int *nbite) {
  double *Ax = (double*)malloc((size_t)n * sizeof(double));
  double *r  = (double*)malloc((size_t)n * sizeof(double));

  double nb = l2_norm(n, b);
  if (nb == 0.0) nb = 1.0;

  int it = 0;
  for (; it < maxit; ++it) {
    dcscmv(n, colptr, rowind, val, x, Ax);
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ax[i];

    double nr = l2_norm(n, r) / nb;
    resvec[it] = nr;
    if (nr < tol) { it++; break; }

    for (int i = 0; i < n; ++i) x[i] += alpha * r[i];
  }

  *nbite = it;
  free(Ax);
  free(r);
}
