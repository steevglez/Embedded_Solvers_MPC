#ifndef QPOASESSW_H
#define QPOASESSW_H

#include <qpOASES.hpp>

void qpoasesSW (elem H[N_QP *N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP]);

#endif // QPOASESSW_H