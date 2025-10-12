#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

typedef struct {
    unsigned int n;     
    unsigned int k;        
    int *offsets;    
    real_t **A;         
    real_t *b;
} MatrizKDiag;

void criaKDiagonal(int n, int k, double **A, double **B);

void genSimetricaPositiva(double *A, double *b, int n, int k, double **ASP, double *bsp, double *tempo);
void geraDLU (double *A, int n, int k, double **D, double **L, double **U, double *tempo);
void geraPreCond(double *D, double *L, double *U, double w, int n, int k, double **M, double *tempo);
double calcResiduoSL (double *A, double *b, double *X, int n, int k, double *tempo);

static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k );
static inline real_t generateRandomB( unsigned int k );

MatrizKDiag* alocaMatrizKDiag(unsigned int n, unsigned int k);
real_t* alocaB(unsigned int n, unsigned int k);

void fillKDiagMatrix(MatrizKDiag *A);
MatrizKDiag* geraMatrizTransposta(const MatrizKDiag *A);
MatrizKDiag *multiplicaMatriz(const MatrizKDiag *A, const MatrizKDiag *B);


#endif // __SISLIN_H__

