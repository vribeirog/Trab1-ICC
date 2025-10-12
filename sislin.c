#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k );
static inline real_t generateRandomB( unsigned int k );

/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return ( (i==j) ? (real_t)(k<<1) : 1.0 )  * (real_t)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomB( unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k<<2) * (real_t)random() * invRandMax;
}

MatrizKDiag* alocaMatrizKDiag(unsigned int n, unsigned int k) {
    MatrizKDiag *A = malloc(sizeof(MatrizKDiag));
    A->n = n;
    A->k = k;

    A->offsets = malloc(k * sizeof(int));
    A->A = malloc(n * sizeof(real_t*));

    int half_band = k / 2;
    for (unsigned int i = 0; i < k; i++) {
        A->offsets[i] = (int)i - half_band;  // e.g., for k=5 → [-2, -1, 0, 1, 2]
    }

    for (unsigned int i = 0; i < n; i++) {
        A->A[i] = calloc(n, sizeof(real_t));
    }

    return A;
}

real_t* alocaB(unsigned int n, unsigned int k) {

  real_t *b = malloc(n * sizeof(real_t));
    
  for (unsigned int i = 0; i < n; i++)
        b[i] = generateRandomB(k);

    return b;
}

void fillKDiagMatrix(MatrizKDiag *M) {
    unsigned int n = M->n;

    // Allocate A and fill with zeros
    M->A = malloc(n * sizeof(real_t *));
    for (unsigned int i = 0; i < n; i++) {
        M->A[i] = calloc(n, sizeof(real_t));
    }

    // Fill each k-diagonal according to offsets
    for (unsigned int d = 0; d < M->k; d++) {
        int offset = M->offsets[d];

        // Positive offsets: upper diagonals
        // Negative offsets: lower diagonals
        for (unsigned int i = 0; i < n; i++) {
            int j = (int)i + offset;
            if (j < 0 || j >= (int)n)
                continue;
            M->A[i][j] = generateRandomA(i, j, M->k);
        }
    }

    // Allocate b and fill with zeros (or random if desired)
    M->b = calloc(n, sizeof(real_t));
}


MatrizKDiag* geraMatrizTransposta(const MatrizKDiag *A) {
    MatrizKDiag *AT = alocaMatrizKDiag(A->n, A->k);
    int k = A->k;

    // da pra melhorar
    for (unsigned int i = 0; i < A->n; i++) {
        for (unsigned int j = 0; j < A->n; j++) {
            AT->A[i][j] = A->A[j][i];
        }
    }

    return AT;
}

MatrizKDiag *multiplicaMatriz(const MatrizKDiag *A, const MatrizKDiag *B) {
    if (A->n != B->n) {
        fprintf(stderr, "Erro: matrizes nao possuem dimensoes compativeis.\n");
        return NULL;
    }

    unsigned int n = A->n;

    // aloca matriz resultado C
    MatrizKDiag *C = malloc(sizeof(MatrizKDiag));
    if (!C) return NULL;

    C->n = n;
    C->k = A->k + B->k - 1;
    C->offsets = NULL;

    // aloca matriz de coefs em C
    C->A = malloc(n * sizeof(real_t *));
    for (unsigned int i = 0; i < n; i++) {
        C->A[i] = calloc(n, sizeof(real_t)); // initialize with 0
    }

    // preencher depois
    C->b = calloc(n, sizeof(real_t));

    // multiplicação de matriz nxn
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            real_t sum = 0.0;
            for (unsigned int k = 0; k < n; k++) {
                sum += A->A[i][k] * B->A[k][j];
            }
            C->A[i][j] = sum;
        }
    }

    return C;
}


/* Cria matriz 'A' k-diagonal e Termos independentes B */
void criaKDiagonal(int n, int k, real_t **A, real_t **B) {

    
  
}

/* Gera matriz simetrica positiva */
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, 
			  real_t **ASP, real_t *bsp, rtime_t *tempo)
{
  *tempo = timestamp();

  *tempo = timestamp() - *tempo;
 
}


void geraDLU (real_t *A, int n, int k,
	      real_t **D, real_t **L, real_t **U, rtime_t *tempo)
{
  *tempo = timestamp();


  *tempo = timestamp() - *tempo;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k,
		 real_t **M, rtime_t *tempo)
{
  *tempo = timestamp();


  *tempo = timestamp() - *tempo;
}


real_t calcResiduoSL (real_t *A, real_t *b, real_t *X,
		      int n, int k, rtime_t *tempo)
{
  *tempo = timestamp();

  real_t *r = calloc(n, sizeof(real_t));

  

  *tempo = timestamp() - *tempo;
}