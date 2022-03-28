#ifndef PDIPSW_H
#define PDIPSW_H
#include <iostream>
#include <cblas.h>

#include "specs.h"
#include "utils.h"
#include "linearSolversSW.h"




#include "mmultSW.h"
using namespace std;
#define MINRES_SW  1
#define CGRAD_SW   2
#define CHOL_SW  3
#define LAPACKCHOL_SW   4
#define LAPACKGESV_SW   5


void pdipSW (elem H[N_QP *N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP], int iterQp, int iterLs, elem tol, int solver);

#endif // PDIPSW_H
