#ifndef PDIPHW_H
#define PDIPHW_H
#include <iostream>

#include "specs.h"
#include "linearSolversHW.h"
#include "mmultHW.h"

using namespace std;


void pdipHW(elem H[N_QP][N_QP], elem h[N_QP], elem M[I_QP][N_QP], elem c[I_QP], elem tk[N_QP], int iterPDIP, int iterLs, elem tol);
#endif
