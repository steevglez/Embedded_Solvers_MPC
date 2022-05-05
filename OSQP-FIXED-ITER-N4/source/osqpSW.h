#ifndef OSQPSW_H
#define OSQPSW_H

#include <osqp.h>
#include "specs.h"
#include "tictoc.h"

void osqpSW(elem H[N_QP *N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP], OSQPWorkspace * & work, OSQPSettings *settings, OSQPData *data,  int *started, tictoc *SOLVERTimer);

#endif // OSQPSW_H
