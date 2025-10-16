#ifndef __GRADCONJ_H__
#define __GRADCONJ_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "matvet.h"
#include "sislin.h"

// Método numérico
real_t gradientesConjugados(real_t **A, real_t *b, real_t *x, int n, real_t tol, int maxit, int* it);

// Função adicional para impressão de resultados
void imprimeResultados(int n, real_t *x, real_t norma, real_t residuo, rtime_t tempo_pc, rtime_t tempo_iter, rtime_t tempo_residuo);

#endif // __GRADCONJ_H__