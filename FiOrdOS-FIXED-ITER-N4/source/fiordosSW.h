#ifndef FIORDOSSW_H
#define FIORDOSSW_H

#include <stdio.h>
#include "specs.h"
#include "tictoc.h"
#include "fiordos4x4_Cparam_solver.h"

void fiordosSW(elem H[N_QP *N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP], fiordos4x4_Cparam_Params *params, fiordos4x4_Cparam_Settings *settings, fiordos4x4_Cparam_Result *result, fiordos4x4_Cparam_Work *work, int *started, tictoc *SOLVERTimer);

#endif // FIORDOSSW_H
