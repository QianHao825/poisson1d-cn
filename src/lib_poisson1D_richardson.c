/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
  int n = *la;

  for (int j = 1; j <= n; ++j) {
    double theta = (double)j * M_PI / (double)(n + 1);
    eigval[j-1] = (2.0 - 2.0 * cos(theta));
  }
}

double eigmax_poisson1D(int *la){
  int n = *la;

  /* j=n gives the maximum */
  double theta = (double)n * M_PI / (double)(n + 1);
  return (2.0 - 2.0 * cos(theta)) ;
}

double eigmin_poisson1D(int *la){
  int n = *la;

  double theta = 1.0 * M_PI / (n + 1.0);   // j = 1
  double lmin = (2.0 - 2.0 * cos(theta)) ;
  return lmin;
}

double richardson_alpha_opt(int *la){
  double lmin = eigmin_poisson1D(la);
  double lmax = eigmax_poisson1D(la);
  return 2.0 / (lmin + lmax);
}

/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  const int n = *la;
  const int ldab = *lab;

  double *Ax = (double*)malloc((size_t)n * sizeof(double));
  double *r  = (double*)malloc((size_t)n * sizeof(double));
  if(!Ax || !r){
    perror("malloc");
    free(Ax); free(r);
    *nbite = 0;
    return;
  }

  double bnorm = cblas_dnrm2(n, RHS, 1);
  if (bnorm == 0.0) bnorm = 1.0;

  int k;
  for(k=0; k<*maxit; ++k){
    /* Ax = A*X */
    cblas_dgbmv(CblasColMajor, CblasNoTrans,n , n, *kl, *ku, 1.0, AB, ldab, X, 1, 0.0, Ax, 1);

    /* r = RHS - Ax */
    cblas_dcopy(n, RHS, 1, r, 1);
    cblas_daxpy(n, -1.0, Ax, 1, r, 1);

    double rnorm = cblas_dnrm2(n, r, 1);
    resvec[k] = rnorm / bnorm;

    if(resvec[k] < *tol){
      ++k; /* number of performed iterations = k */
      break;
    }

    /* X = X + alpha*r */
    cblas_daxpy(n, *alpha_rich, r, 1, X, 1);
  }

  *nbite = (k > *maxit) ? *maxit : k;

  free(Ax);
  free(r);
}

/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){

  (void)kv; /* kv unused here: we assume standard GB for mat-vec */
  const int n = *la;
  const int ldab = *lab;
  const int diag = *ku;

  /* zero fill */
  for(int j=0; j<n; ++j){
    for(int i=0; i<ldab; ++i){
      MB[indexABCol(i,j,lab)] = 0.0;
    }
  }

  for(int j=0; j<n; ++j){
    MB[indexABCol(diag, j, lab)] = AB[indexABCol(diag, j, lab)];
  }

}

/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  (void)kv;
  const int n = *la;
  const int ldab = *lab;
  const int diag = *ku;
  const int sub  = diag + 1; /* one sub diagonal for tridiag */

  /* zero fill */
  for(int j=0; j<n; ++j){
    for(int i=0; i<ldab; ++i){
      MB[indexABCol(i,j,lab)] = 0.0;
    }
  }

  for(int j=0; j<n; ++j){
    /* diag */
    MB[indexABCol(diag, j, lab)] = AB[indexABCol(diag, j, lab)];
    /* sub (exists for columns 0..n-2 in GB storage) */
    if(sub < ldab){
      MB[indexABCol(sub, j, lab)] = AB[indexABCol(sub, j, lab)];
    }
  }
}

/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  const int n = *la;
  const int ldab = *lab;

  double *Ax = (double*)malloc((size_t)n * sizeof(double));
  double *r  = (double*)malloc((size_t)n * sizeof(double));
  double *z  = (double*)malloc((size_t)n * sizeof(double));
  if(!Ax || !r || !z){
    perror("malloc");
    free(Ax); free(r); free(z);
    *nbite = 0;
    return;
  }

  double bnorm = cblas_dnrm2(n, RHS, 1);
  if (bnorm == 0.0) bnorm = 1.0;

  /* detect whether MB is diagonal-only (Jacobi) or lower+diag (GS) */
  const int diag = *ku;
  const int sub  = diag + 1;

  int mb_has_sub = 0;
  if(sub < ldab){
    for(int j=0; j<n; ++j){
      if(MB[indexABCol(sub, j, lab)] != 0.0){
        mb_has_sub = 1;
        break;
      }
    }
  }

  int k;
  for(k=0; k<*maxit; ++k){
    /* Ax = A*X */
    cblas_dgbmv(CblasColMajor, CblasNoTrans, n, n, *kl, *ku, 1.0, AB, ldab, X, 1, 0.0, Ax, 1);

    /* r = RHS - Ax */
    cblas_dcopy(n, RHS, 1, r, 1);
    cblas_daxpy(n, -1.0, Ax, 1, r, 1);

    double rnorm = cblas_dnrm2(n, r, 1);
    resvec[k] = rnorm / bnorm;

    if(resvec[k] < *tol){
      ++k;
      break;
    }

    /* Solve MB * z = r */
    if(!mb_has_sub){
      /* Jacobi: z_i = r_i / D_i */
      for(int i=0; i<n; ++i){
        double d = MB[indexABCol(diag, i, lab)];
        z[i] = (d != 0.0) ? (r[i] / d) : r[i];
      }
    }else{
      /* Gauss-Seidel: forward substitution on lower bidiagonal */
      for(int i=0; i<n; ++i){
        double d = MB[indexABCol(diag, i, lab)];
        double rhs = r[i];
        if(i > 0){
          double l = MB[indexABCol(sub, i-1, lab)]; /* A_{i,i-1} stored in column i-1 */
          rhs -= l * z[i-1];
        }
        z[i] = (d != 0.0) ? (rhs / d) : rhs;
      }
    }

    /* X = X + z */
    cblas_daxpy(n, 1.0, z, 1, X, 1);
  }

  *nbite = (k > *maxit) ? *maxit : k;

  free(Ax);
  free(r);
  free(z);
}
