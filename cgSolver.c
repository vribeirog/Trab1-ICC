// Isadora Botassari - GRR20206872
// Victor Ribeiro Garcia - GRR20203954

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "utils.h"
#include "matvet.h"
#include "sislin.h"
#include "gradconj.h"

int main() {
    srandom(20252);

    int n, k, maxit;
    real_t omega, epsilon, norma = 0.0, residuo = 0.0;
    real_t **A, *b, *x, **ASP, *bsp, **D, **L, **U, **M;
    rtime_t tempo_simetrica = 0.0, tempo_dlu = 0.0, tempo_pc_parcial = 0.0, tempo_iter = 0.0, tempo_residuo = 0.0;

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

    // Alocar todos os vetores necessários
    if (aloca_vetores(&b, &bsp, &x, n) != 0) {
        return 1;
    }
    
    // Alocar todas as matrizes necessárias
    if (aloca_matrizes(&A, &ASP, &D, &L, &U, &M, n) != 0) {
        free_all(&A, &b, &x, &ASP, &bsp, &D, &L, &U, &M);
        return 1;
    }

    // Gerar matriz A k-diagonal e vetor de termos independetes b
    criaKDiagonal(n, k, A, b);

    // Gerar matriz simétrica positiva usando a função genSimetricaPositiva
    // ASP = A * AT
    // bsp = AT * b
    genSimetricaPositiva(A, b, n, k, ASP, bsp, &tempo_simetrica);

    // Gerar matrizes D, L e U para o pré-condicionador de Jacobi, se necessário
    if (omega >= 0.0)
        geraDLU(ASP, n, k, D, L, U, &tempo_dlu); 

    // Gerar matriz M⁻¹ de acordo com o pré-condicionador
    geraPreCond(D, L, U, omega, n, k, M, &tempo_pc_parcial); 
    
    // Resolver o sistema linear usando gradientes conjugados (com ou sem pré-condicionador)
    if (omega == 0.0)
        norma = gradientesConjugadosPrecond(M, ASP, bsp, x, n, epsilon, maxit, &tempo_iter);
    else if (omega == -1.0)
        norma = gradientesConjugados(ASP, bsp, x, n, epsilon, maxit, &tempo_iter);
    else {
        fprintf(stderr, "Erro: valor de omega inválido. Use -1.0 (sem pré-condicionador) ou 0.0 (Jacobi).\n");
        free_all(&A, &b, &x, &ASP, &bsp, &D, &L, &U, &M);
        return 1;
    }

    residuo = calcResiduoSL(ASP, bsp, x, n, k, &tempo_residuo);
    
    imprimeResultados(n, x, norma, residuo, tempo_simetrica + tempo_dlu + tempo_pc_parcial, tempo_iter, tempo_residuo);
    

    // Liberar toda a memoria alocada
    free_all(&A, &b, &x, &ASP, &bsp, &D, &L, &U, &M);

    return 0;
}