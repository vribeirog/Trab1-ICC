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

/* Cria matriz 'A' k-diagonal e Termos independentes B */
void criaKDiagonal(int n, int k, real_t **A, real_t *B) {
    int banda = k / 2; // Número de diagonais acima (e abaixo) da diagonal principal

    // Inicializar matriz A com zeros
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = 0.0;
        }
    }

    // Preencher as k diagonais da matriz
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Verificar se o elemento está dentro da banda k-diagonal
            if (abs(i - j) <= banda) {
                A[i][j] = generateRandomA(i, j, k);
            }
        }
    }

    // Preencher vetor B
    for (int i = 0; i < n; i++) {
        B[i] = generateRandomB(k);
    }
}

/* Gera matriz simetrica positiva */
void genSimetricaPositiva(real_t **A, real_t *b, int n, int k, 
			  real_t **ASP, real_t *bsp, rtime_t *tempo)
{
    *tempo = timestamp();
    
    // Alocar matriz transposta AT usando a mesma estrutura contígua
    real_t **AT = aloca_matriz(n, 0);
    if (!AT) {
        fprintf(stderr, "Erro ao alocar matriz transposta AT na genSimetricaPositiva.\n");
        exit(1);
    }
    
    // Calcular transposta: AT[j][i] = A[i][j]
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            AT[j][i] = A[i][j];
        }
    }
    
    printf("matriz transposta:\n");
    imprime_matriz(AT, n);
    printf("\n");

    // A matriz ASP já foi alocada pela função aloca_matrizes
    // Calcular ASP = A * AT (matriz simétrica positiva definida)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            real_t sum = 0.0;
            for (int k_idx = 0; k_idx < n; k_idx++) {
                sum += A[i][k_idx] * AT[k_idx][j];
            }
            ASP[i][j] = sum;
        }
    }
    
    // Calcular bsp = AT * b (lado direito do sistema simétrico)
    for (int i = 0; i < n; i++) {
        bsp[i] = 0.0;
        for (int j = 0; j < n; j++) {
            bsp[i] += AT[i][j] * b[j];
        }
    }
    
    // Liberar matriz transposta temporária
    free_matriz(AT);
    
    *tempo = timestamp() - *tempo;
}

// Preencher D, L, U
// D: diagonal principal de A
// L: parte inferior de A (sem diagonal principal)
// U: parte superior de A (sem diagonal principal)
void geraDLU (real_t **A, int n, int k,
	      real_t **D, real_t **L, real_t **U, rtime_t *tempo)
{
    *tempo = timestamp();

    int banda = k / 2; // Número de diagonais acima (e abaixo) da diagonal principal
    int i, j;

    // U e L devem ter a diagonal principal nula.
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j <= banda + i && j < n; j++)
            U[i][j] = A[i][j];

    for (int i = 0; i < n; i++)
        for (int j = i - 1; j >= i - banda && j >= 0; j--)
            L[i][j] = A[i][j];
    
    for (i = 0; i < n; i++)
        D[i][i] = A[i][i];

    *tempo = timestamp() - *tempo;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(real_t **D, real_t **L, real_t **U, real_t w, int n, int k,
		 real_t **M, rtime_t *tempo)
{
    *tempo = timestamp();

    // Se w == -1.0, sem pré-condicionador (M = I)
    if (w == -1.0){
         for (int i = 0; i < n; i++)
            M[i][i] = 1.0;
    }
    else if (w == 0.0){
        // Pré-condicionador de Jacobi (M = D)
        for (int i = 0; i < n; i++)
            M[i][i] = D[i][i];
    }

    *tempo = timestamp() - *tempo;
}


real_t calcResiduoSL (real_t **A, real_t *b, real_t *X,
		      int n, int k, rtime_t *tempo)
{
  *tempo = timestamp();

  *tempo = timestamp() - *tempo;
}

// ==================== FUNÇÕES MATRIZES E VETORES ====================

// Imprime matriz
void imprime_matriz(real_t **mat, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%6.3f ", mat[i][j]);
        }
        printf("\n");
    }
}

// Imprime vetor
void imprime_vetor(real_t *vec, int n) {
    for (int i = 0; i < n; i++) {
        printf("%6.3f ", vec[i]);
    }
    printf("\n");
}

// Aloca uma matriz bidimensional de forma contígua na memória
// Permite acesso mat[i][j] com apenas uma alocação
// Se zero_init == 1 inicializa todos os elementos com zero
real_t** aloca_matriz(int n, int zero_init) {
    real_t **mat;
    
    // Aloca um vetor com os ponteiros e os elementos da matriz
    if (zero_init) {
        mat = calloc(n * sizeof(real_t*) + n * n * sizeof(real_t), 1);
    } else {
        mat = malloc(n * sizeof(real_t*) + n * n * sizeof(real_t));
    }
    
    if (!mat) return NULL;
    
    // Ajusta o ponteiro da primeira linha
    mat[0] = (real_t*)(mat + n);
    
    // Ajusta os ponteiros das demais linhas (i > 0)
    for (int i = 1; i < n; i++) {
        mat[i] = mat[0] + (i * n);
    }
    
    return mat;
}

// Libera uma matriz alocada com aloca_matriz
void free_matriz(real_t **mat) {
    if (mat) {
        free(mat);
    }
}

// Aloca todos os vetores necessários para o programa
int aloca_vetores(real_t **b, real_t **bsp, real_t **x, int n) {
    // Alocar vetor b
    *b = malloc(n * sizeof(real_t));
    if (!*b) {
        fprintf(stderr, "Erro ao alocar memoria para vetor b.\n");
        return 1;
    }
    
    // Alocar vetor bsp (AT * b)
    *bsp = malloc(n * sizeof(real_t));
    if (!*bsp) {
        fprintf(stderr, "Erro ao alocar memoria para vetor bsp.\n");
        free(*b);
        return 1;
    }
    
    // Alocar vetor solução x
    *x = calloc(n, sizeof(real_t)); // Inicializado com zeros
    if (!*x) {
        fprintf(stderr, "Erro ao alocar memoria para vetor solução x.\n");
        free(*b);
        free(*bsp);
        return 1;
    }
    
    return 0; // Sucesso
}

// Aloca todas as matrizes necessárias para o programa
int aloca_matrizes(real_t ***A, real_t ***ASP, real_t ***D, real_t ***L, real_t ***U, real_t ***M, int n) {
    // Alocar matriz A (k-diagonal original)
    *A = aloca_matriz(n, 0);
    if (!*A) {
        fprintf(stderr, "Erro ao alocar memoria para matriz A.\n");
        return 1;
    }
    
    // Alocar matriz ASP (simétrica positiva)
    *ASP = aloca_matriz(n, 0);
    if (!*ASP) {
        fprintf(stderr, "Erro ao alocar memoria para matriz ASP.\n");
        free_matriz(*A);
        return 1;
    }
    
    // Alocar matriz D (diagonal)
    *D = aloca_matriz(n, 1);
    if (!*D) {
        fprintf(stderr, "Erro ao alocar memoria para matriz D.\n");
        free_matriz(*A);
        free_matriz(*ASP);
        return 1;
    }
    
    // Alocar matriz L (triangular inferior)
    *L = aloca_matriz(n, 1);
    if (!*L) {
        fprintf(stderr, "Erro ao alocar memoria para matriz L.\n");
        free_matriz(*A);
        free_matriz(*ASP);
        free_matriz(*D);
        return 1;
    }
    
    // Alocar matriz U (triangular superior)
    *U = aloca_matriz(n, 1);
    if (!*U) {
        fprintf(stderr, "Erro ao alocar memoria para matriz U.\n");
        free_matriz(*A);
        free_matriz(*ASP);
        free_matriz(*D);
        free_matriz(*L);
        return 1;
    }
    
    // Alocar matriz M (pré-condicionadora)
    *M = aloca_matriz(n, 1);
    if (!*M) {
        fprintf(stderr, "Erro ao alocar memoria para matriz pré-condicionadora M.\n");
        free_matriz(*A);
        free_matriz(*ASP);
        free_matriz(*D);
        free_matriz(*L);
        free_matriz(*U);
        return 1;
    }

    return 0;
}

// Libera toda a memória alocada para vetores e matrizes
void free_all(real_t ***A, real_t **b, real_t **x, real_t ***ASP, real_t **bsp, 
              real_t ***D, real_t ***L, real_t ***U, real_t ***M) {
    
    // Liberar vetores
    if (b && *b) {
        free(*b);
        *b = NULL;
    }
    
    if (bsp && *bsp) {
        free(*bsp);
        *bsp = NULL;
    }
    
    if (x && *x) {
        free(*x);
        *x = NULL;
    }
    
    // Liberar matrizes bidimensionais
    if (A && *A) {
        free_matriz(*A);
        *A = NULL;
    }
    
    if (ASP && *ASP) {
        free_matriz(*ASP);
        *ASP = NULL;
    }
    
    if (D && *D) {
        free_matriz(*D);
        *D = NULL;
    }
    
    if (L && *L) {
        free_matriz(*L);
        *L = NULL;
    }
    
    if (U && *U) {
        free_matriz(*U);
        *U = NULL;
    }
    
    if (M && *M) {
        free_matriz(*M);
        *M = NULL;
    }
}