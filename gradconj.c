#include <stdio.h>
#include <stdlib.h>    
#include <string.h>
#include <math.h>

#include "utils.h"
#include "matvet.h"
#include "sislin.h"
#include "gradconj.h"

// Calcula produtor escalar (dot product) entre vetores
real_t dot(real_t *a, real_t *b, int n) {
    real_t s = 0.0;
    for (int i = 0; i < n; i++) {
        s += a[i] * b[i];
    }
    return s;
}

// Calcula norma de um vetor a
real_t norma(real_t *a, int n) {
    return sqrt(dot(a, a, n));
}

// Calcula produto entre uma matriz e um vetor de dimensoes compativeis
// talvez add erro caso dimensoes nao sejam compativeis
void prodMatVet(real_t **A, real_t *x, real_t *y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) {
            y[i] += A[i][j] * x[j];
        }
    }
}

// Método numérico de gradientes conjugados
real_t gradientesConjugados(real_t **A, real_t *b, real_t *x, int n, real_t tol, int maxit, int* it) {
    real_t *residuo = malloc(n * sizeof(real_t));
    if (!residuo) {
        fprintf(stderr, "Erro ao alocar vetor residuo na gradientesConjugados.\n");
        exit(1);
    }

    real_t *search_direction = malloc(n * sizeof(real_t));
    if (!search_direction) {
        fprintf(stderr, "Erro ao alocar vetor search_direction na gradientesConjugados.\n");
        free(residuo);
        exit(1);
    }

    real_t *A_search_direction = malloc(n * sizeof(real_t));
    if (!A_search_direction) {
        fprintf(stderr, "Erro ao alocar vetor A_search_direction na gradientesConjugados.\n");
        free(residuo);
        free(search_direction);
        exit(1);
    }

    int iter = 0;

    // calculo do residuo = b - A*x
    prodMatVet(A, x, residuo, n);
    for (int i = 0; i < n; i++) {
        residuo[i] = b[i] - residuo[i];
        search_direction[i] = residuo[i];
    }

    real_t old_resid_norm = norma(residuo, n);

    // itera enquanto a norma é menor que a tolerancia (epsilon) ou nao ultrapassa maxit
    while ((old_resid_norm > tol) && (iter < maxit)) {
        prodMatVet(A, search_direction, A_search_direction, n);
        real_t denom = dot(search_direction, A_search_direction, n);
        real_t step_size = (old_resid_norm * old_resid_norm) / denom;

        // x = x + step_size * search_direction
        for (int i = 0; i < n; i++) {
            x[i] += step_size * search_direction[i];
            residuo[i] -= step_size * A_search_direction[i];
        }

        real_t new_resid_norm = norma(residuo, n);
        real_t beta = (new_resid_norm * new_resid_norm) / (old_resid_norm * old_resid_norm);

        for (int i = 0; i < n; i++) {
            search_direction[i] = residuo[i] + beta * search_direction[i];
        }

        old_resid_norm = new_resid_norm;
        iter++;
    }

    *it = iter;

    free(residuo);
    free(search_direction);
    free(A_search_direction);
}


// Imprime resultados onde:
// - n: tamanho do vetor solução x;
// - x_1 x_2 ... x_n: valores do vetor solução x;
// - norma: norma máxima do erro aproximado em x após última iteração;
// - resíduo: norma euclidiana do resíduo;
// - tempo_pc: tempo para calcular a matriz pré-condicionante M e preparar o SL para o uso do pré-condicionante;
// - tempo_iter: tempo médio para calcular uma iteração do método, inclusive o cálculo do erro;
// - tempo_residuo: tempo para calcular a norma euclidiana do resíduo ao final do processo.
void imprimeResultados(int n, real_t *x, real_t norma, real_t residuo, rtime_t tempo_pc, rtime_t tempo_iter, rtime_t tempo_residuo){
    printf("%d\n", n);

    if (!x){
        printf("Vetor solucoes x é nulo!\n");
        return;
    }
    for(int i=0; i < n-1; i++){
        printf("%.16g ", x[i]);
    }
    printf("%.16g\n", x[n-1]);

    printf("%.8g\n", norma);
    printf("%.16g\n", residuo);
    printf("%.8g\n", tempo_pc);
    printf("%.8g\n", tempo_iter);
    printf("%.8g\n", tempo_residuo);
}