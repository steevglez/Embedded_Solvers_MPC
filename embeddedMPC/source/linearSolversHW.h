#ifndef LINEARSOLVERSHW_H
#define LINEARSOLVERSHW_H

#include "specs.h"
#include <stdio.h>
#include <math.h>


// linear solvers
# pragma SDS data zero_copy(Ak[0:N_QP*N_QP])
void minresHW(elem Ak[N_QP][N_QP], elem bk[N_QP], elem zko[N_QP], elem zk[N_QP], elem iter, elem tol);

void cgradHW(elem Ak[N_QP][N_QP], elem bk[N_QP], elem zko[N_QP], elem zk[N_QP], elem tol);

void cholHW(elem Ak[N_QP][N_QP], elem bk[N_QP], elem zk[N_QP]);

#endif
