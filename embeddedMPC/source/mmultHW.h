#ifndef MMULTHW_H
#define MMULTHW_H

#include "specs.h"
#define MAX_SIZE I_QP



elem adderTree (elem input[MAX_SIZE]);
#pragma SDS data zero_copy(A[0:m*n], B[0:n*p], C[0:m*n])
#pragma SDS data access_pattern(A:SEQUENTIAL, B:SEQUENTIAL, C:SEQUENTIAL)
void mmultHW (elem A[MAX_SIZE][MAX_SIZE], elem B[MAX_SIZE][MAX_SIZE], elem C[MAX_SIZE][MAX_SIZE], int m, int n, int p);
void mmult_N_I_N (elem A[N_QP][I_QP], elem B[I_QP][N_QP], elem C[N_QP][N_QP]);
void mmult_N_N_1 (elem A[N_QP][N_QP], elem B[N_QP], elem C[N_QP]);
void mmult_N_I_1 (elem A[N_QP][I_QP], elem B[I_QP], elem C[N_QP]);
void mmult_I_N_1 (elem A[I_QP][N_QP], elem B[N_QP], elem C[I_QP]);

#endif
