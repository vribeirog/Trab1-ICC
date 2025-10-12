#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "utils.h"
#include "sislin.h"

int main() {
    srandom(20252);

    int n, k, maxit;
    real_t omega, epsilon;

    scanf("%d", &n);
    scanf("%d", &k);
    scanf("%lf", &omega);
    scanf("%d", &maxit);
    scanf("%lf", &epsilon);

    MatrizKDiag *A = alocaMatrizKDiag(n, k) ;
    
    fillKDiagMatrix(A);
    printf("matriz orig:\n");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) 
            printf("%6.3f ", A->A[i][j]);
        printf("\n");
    }
    
    printf("\n");

    real_t *b = alocaB(n,k);

    MatrizKDiag *AT = geraMatrizTransposta(A);

    printf("matriz transposta:\n");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) 
            printf("%6.3f ", AT->A[i][j]);
        printf("\n");
    }

    MatrizKDiag *C = multiplicaMatriz(A, AT);
     printf("\n");
    
    printf("matriz simetrica:\n");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) 
            printf("%6.3f ", C->A[i][j]);
        printf("\n");
    }


    return 0;
}