#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

void criaKDiagonal(int n, int k, real_t **A, real_t *B);

void genSimetricaPositiva(real_t **A, real_t *b, int n, int k, real_t **ASP, real_t *bsp, rtime_t *tempo);
void geraDLU (real_t **A, int n, int k, real_t **D, real_t **L, real_t **U, rtime_t *tempo);
void geraPreCond(real_t **D, real_t **L, real_t **U, real_t w, int n, int k, real_t **M, rtime_t *tempo);
real_t calcResiduoSL (real_t **A, real_t *b, real_t *X, int n, int k, rtime_t *tempo);

static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k );
static inline real_t generateRandomB( unsigned int k );

// Impressão
void imprime_matriz(real_t **mat, int n);
void imprime_vetor(real_t *vet, int n);

// Funções de alocação e desalocação
int aloca_vetores(real_t **b, real_t **bsp, real_t **x, int n);
int aloca_matrizes(real_t ***A, real_t ***ASP, real_t ***D, real_t ***L, real_t ***U, real_t ***M, int n);
void free_all(real_t ***A, real_t **b, real_t **x, real_t ***ASP, real_t **bsp, real_t ***D, real_t ***L, real_t ***U, real_t ***M);

// Funções auxiliares para alocação/liberação de matrizes
real_t** aloca_matriz(int n, int zero_init);
void free_matriz(real_t **mat);

// Método numérico
void gradientesConjugados(real_t **A, real_t *b, real_t *x, int n, real_t tol, int maxit, int* it);

#endif // __SISLIN_H__

