#include <stdio.h>
#include <stdlib.h>
#include "lib_poisson1D.h"

int main() {
    int n   = 6;  
    int la  = n;
    int kv  = 1;
    int kl  = 1;
    int ku  = 1;
    int lab = kv + kl + ku + 1;   // = 4

    double *AB = malloc(lab * la * sizeof(double));

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

    printf("AB matrix (n rows, 4 columns):\n\n");

    for (int i = 0; i < la; i++) {      
        for (int j = 0; j < lab; j++) {  
            printf("%8.3f  ", AB[i * lab + j]);
        }
        printf("\n");
    }

    free(AB);
    return 0;
}
