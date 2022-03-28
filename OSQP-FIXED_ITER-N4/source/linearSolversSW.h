#ifndef LINEARSOLVERSSW_H
#define LINEARSOLVERSSW_H
#include <math.h>
#include <cblas.h>
#include <iostream>
#include "specs.h"
#include "utils.h"

// https://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=4221
#include <complex.h>
#define lapack_complex_float    float _Complex
#define lapack_complex_double   double _Complex
#include <lapacke.h>
void minresSW(elem* Ak, elem* bk, elem* zko, elem* zk, elem iter, elem tol);

void cgradSW(elem* Ak, elem* bk, elem* zko, elem* zk, elem tol);

void cholSW(elem* Ak, elem* bk, elem* zk);

void lapackcholSW(elem* Ak, elem* bk, elem* zk);

void lapackgesvSW(elem* Ak, elem* bk, elem* zk);

#endif
