#ifndef OSQPSW_H
#define OSQPSW_H

#include "osqp.h"

using namespace std;

void osqpSW(elem H[N_QP *N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP]);

#endif // OSQPSW_H