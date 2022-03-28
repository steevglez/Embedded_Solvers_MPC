#ifndef MPCHW_H
#define MPCHW_H

#include<iostream>

#include"specs.h"
#include"pdipHW.h"

using namespace std;

void mpcHW(elem x[N_SYS], elem yref[P_SYS], elem xmin[N_SYS], elem xmax[N_SYS], elem umin[M_SYS], elem umax[M_SYS], elem Acal[N_SYS*N_QP][N_SYS], elem AcalQOcal[N_SYS][N_QP], elem H[N_QP][N_QP], elem Mx[N_SYS*N_QP*2+2*N_QP][N_QP], elem L_col[N_SYS+P_SYS], int iterPDIP, int iterLS, elem u[N_SYS]);

void get_h(elem AcalQOcal[N_SYS][N_QP], elem x_tilde[N_SYS], elem h[N_QP]);
void get_cx(elem Acal[N_SYS*N_QP][N_SYS], elem xmax[N_SYS], elem xmin[N_SYS], elem x_star[N_SYS], elem x_tilde[N_SYS], elem a[N_QP], elem b[N_QP], elem cx[6*N_QP]);

#endif
