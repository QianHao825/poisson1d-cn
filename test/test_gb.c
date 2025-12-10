#include <stdio.h>
#include <stdlib.h>
#include "lib_poisson1D.h"

int main() {
    int n   = 6;     // 内部点个数
    int la  = n;
    int kv  = 1;
    int kl  = 1;
    int ku  = 1;
    int lab = kv + kl + ku + 1;   // = 4

    double *AB = malloc(lab * la * sizeof(double));

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

    printf("AB matrix (n rows, 4 columns):\n\n");

    // 和 write_GB_operator_colMajor_poisson1D 一样：按“行”打印
    for (int i = 0; i < la; i++) {       // i: 0..n-1 → 行
        for (int j = 0; j < lab; j++) {  // j: 0..3   → 列
            printf("%8.3f  ", AB[i * lab + j]);
        }
        printf("\n");
    }

    free(AB);
    return 0;
}
