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
    // ASP = A * AT
    // bsp = AT * b
    real_t **A, *b, *x, **ASP, *bsp, **D, **L, **U, **M;
    rtime_t tempo_simetrica, tempo_dlu, tempo_pc = 0.0;

    scanf("%d", &n); // Dimensao da matriz (n x n)
    scanf("%d", &k); // Numero de diagonais da matriz A
    // k deve ser maior que 1 e impar
    if (k < 1 || k % 2 == 0) {
        fprintf(stderr, "Erro: k deve ser maior que 1 e impar.\n");
        return 1;
    }
    scanf("%lf", &omega); // Pre-condicionador: omega = -1.0 (sem pre-condicionador) e omega = 0.0 (Jacobi)
    scanf("%d", &maxit); // Numero maximo de iteracoes a serem executadas
    scanf("%lf", &epsilon); // Erro aproximado absoluto maximo tolerado

    printf("Parâmetros: omega=%.3f, maxit=%d, epsilon=%e\n", omega, maxit, epsilon);

    // Alocar todos os vetores necessários
    if (aloca_vetores(&b, &bsp, &x, n) != 0) {
        return 1;
    }
    
    // Alocar todas as matrizes necessárias
    if (aloca_matrizes(&A, &ASP, &D, &L, &U, &M, n) != 0) {
        free_all(&A, &b, &x, &ASP, &bsp, &D, &L, &U, &M);
        return 1;
    }

    criaKDiagonal(n, k, A, b);
    
    printf("matriz orig:\n");
    imprime_matriz(A, n);
    printf("\n");
    
    printf("vetor b:\n");
    imprime_vetor(b, n);
    printf("\n");

    // Gerar matriz simétrica positiva usando a função genSimetricaPositiva
    genSimetricaPositiva(A, b, n, k, ASP, bsp, &tempo_simetrica);
    
    printf("matriz simetrica:\n");
    imprime_matriz(ASP, n);
    printf("\n");
    
    printf("vetor AT * b:\n");
    imprime_vetor(bsp, n);
    printf("\n");

    if (omega >= 0.0)
        geraDLU(A, n, k, D, L, U, &tempo_dlu); 

    imprime_matriz(D, n);
    printf("\n");
    imprime_matriz(L, n);
    printf("\n");
    imprime_matriz(U, n);
    printf("\n");

    geraPreCond(D, L, U, omega, n, k, M, &tempo_pc); 
    

    // Liberar toda a memoria alocada
    free_all(&A, &b, &x, &ASP, &bsp, &D, &L, &U, &M);

    return 0;
}