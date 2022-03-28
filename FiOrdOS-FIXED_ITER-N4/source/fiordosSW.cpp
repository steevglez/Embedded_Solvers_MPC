// Compilar: gcc fiordos_Cparam_main.c fiordos_Cparam_solver.c -lm -o fiordos_Cparam
// Ejecutar: ./fiordos_Cparam

#include "fiordosSW.h"


void fiordosSW(elem H[N_QP *N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP], fiordos4x4_Cparam_Params *params, fiordos4x4_Cparam_Settings *settings, fiordos4x4_Cparam_Result *result, fiordos4x4_Cparam_Work *work, int *started, tictoc *SOLVERTimer){

    int i;
    if(*started == 0){
    	/* define all parametric data */
        for(i = 0 ; i < N_QP*N_QP ; i++){
        	params->H[i] = H[i];
    	}
    	for(i = 0 ; i < N_QP ; i++){
      		params->g[i] = h[i];
    	}
    	for(i = 0 ; i < I_QP*N_QP ; i++){
      		params->Ai[i] = M[i];
    	}
    	for(i = 0 ; i < I_QP ; i++){
      		params->bi[i] = c[i];
    	}
    	*started = 1;
    }else if(*started == 1){
    	/* define some parametric data and load previous results */
    	for(i = 0 ; i < N_QP ; i++){
      		params->g[i] = h[i];
      		settings->algoInner.init[i] = result->x[i];
    	}
        for(i = 0 ; i < I_QP ; i++){
      		params->bi[i] = c[i];
      		settings->algoOuter.init[i] = result->la[i];
    	}	
    }
        
    /* solve the problem */
    SOLVERTimer->tic();
    fiordos4x4_Cparam_solve(params, settings, result, work);
    SOLVERTimer->toc();

    for(i = 0 ; i < N_QP ; i++){
      tk[i] = result->x[i];
    }

}
