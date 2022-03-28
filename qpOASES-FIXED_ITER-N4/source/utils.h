#ifndef LS_UTILS
#define LS_UTILS

#include "specs.h"
#include <cblas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// BLAS wrappers
#if defined ELEM_DOUBLE
//#define mmultSW(A,transA, B, transB, C, m, n, p) cblas_dgemm(CblasRowMajor, transA, transB, m, p, n, 1, A, n, B, p, 0, C, p)
#define POSV(N, A, B) LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', N, 1, A, N, B, 1)
#define GESV(N, A, B, ipiv) LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, 1, A, N, ipiv, B, 1)
#define COPY(N, A, B) cblas_dcopy(N, A, 1, B, 1)
#define DOT(n, x, incx, y, incy) cblas_ddot(n, x, incx, y, incy)
#define AXPY(n, a, x, incx, y, incy) cblas_daxpy(n, a, x, incx, y, incy)
#elif defined ELEM_FLOAT
//#define mmultSW(A,transA, B, transB, C, m, n, p) cblas_sgemm(CblasRowMajor, transA, transB, m, p, n, 1, A, n, B, p, 0, C, p)
#define POSV(N, A, B) LAPACKE_sposv(LAPACK_ROW_MAJOR, 'U', N, 1, A, N, B, 1)
#define GESV(N, A, B, ipiv) LAPACKE_sgesv(LAPACK_ROW_MAJOR, N, 1, A, N, ipiv, B, 1)
#define COPY(N, A, B) cblas_scopy(N, A, 1, B, 1)
#define DOT(n, x, incx, y, incy) cblas_sdot(n, x, incx, y, incy)
#define AXPY(n, a, x, incx, y, incy) cblas_saxpy(n, a, x, incx, y, incy)
#else
#error Something is wrong with how elem is defined
#endif

void displayMatrix(double* C, int rows, int cols);
void displayMatrix(elem* C, int rows, int cols);
void genRandArray(double min, double max, int size, double *array);
void genSymPosMat(double max, int size, double *matrix);
int compare(double* gold, elem* result, int size, double th);
int compare(float* gold, elem* result, int size, double th);
void cast2elem(double *inputArray, elem *outputArray, int size);
void cast2float(double *inputArray, float *outputArray, int size);


#endif
