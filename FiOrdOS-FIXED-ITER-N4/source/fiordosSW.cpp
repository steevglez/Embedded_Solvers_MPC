#include "fiordosSW.h"

void fiordosSW(elem H[N_QP *N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP], fiordos4x4_Cparam_Params *params, fiordos4x4_Cparam_Settings *settings, fiordos4x4_Cparam_Result *result, fiordos4x4_Cparam_Work *work, int *started, tictoc *SOLVERTimer){

    int i;
    
    // Si el codigo aun no ha realizado una iteracion, se cargan los todos los datos
    // a las estructuras de datos al solver
    if(*started == 0){
    	// Matriz H - Hessiana
        for(i = 0 ; i < N_QP*N_QP ; i++){
        	params->H[i] = H[i];
    	}
    	// Vector h - Funcional o lineal
    	for(i = 0 ; i < N_QP ; i++){
      		params->g[i] = h[i];
    	}
    	// Matriz M(Mx) - Matriz de restricciones
    	for(i = 0 ; i < I_QP*N_QP ; i++){
      		params->Ai[i] = M[i];
    	}
    	// Vector c - Vector de restricciones
    	for(i = 0 ; i < I_QP ; i++){
      		params->bi[i] = c[i];
    	}
    	*started = 1; 	// Marca que ya se ha hecho una iteracion
    }
    // Si el codigo ya ha hecho una iteracion, al solver solo se cargan los datos 
    // de los vectores h y c
    else if(*started == 1){
    	// Actualizacion de vector h
    	for(i = 0 ; i < N_QP ; i++){
      		params->g[i] = h[i];
      		settings->algoInner.init[i] = result->x[i];
    	}
    	// Actualizacion de vector c
        for(i = 0 ; i < I_QP ; i++){
      		params->bi[i] = c[i];
      		settings->algoOuter.init[i] = result->la[i];
    	}	
    }
        
    // Resolver el problema
    SOLVERTimer->tic();
    fiordos4x4_Cparam_solve(params, settings, result, work);
    SOLVERTimer->toc();

    // Se almacena la solucion obtenida 
    for(i = 0 ; i < N_QP ; i++){
      tk[i] = result->x[i];
    }

}
