#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib_poisson1D.h"

int main(void)
{
    int n = 100;
    int maxit = 1000;
    double tol = 1e-6;
    double alpha = 0.5;

    int nnz = nnz_poisson1D(n);

    int *rowptr = malloc((n+1)*sizeof(int));
    int *colind = malloc(nnz*sizeof(int));
    double *val = malloc(nnz*sizeof(double));

    double *b = malloc(n*sizeof(double));
    double *x = calloc(n, sizeof(double));
    double *resvec = malloc(maxit*sizeof(double));
    int nbite;

    /* Build CSR matrix */
    set_CSR_poisson1D(n, rowptr, colind, val);

    /* RHS f = 1 */
    for (int i = 0; i < n; ++i)
        b[i] = 1.0;

    /* Richardson CSR */
    richardson_alpha_csr(n, rowptr, colind, val,
                          b, x, alpha, tol, maxit,
                          resvec, &nbite);

    printf("EX10 (CSR Richardson): %d iterations\n", nbite);
    printf("last residual = %.6e\n", resvec[nbite-1]);


    free(rowptr);
    free(colind);
    free(val);
    free(b);
    free(x);
    free(resvec);

    return 0;
}
