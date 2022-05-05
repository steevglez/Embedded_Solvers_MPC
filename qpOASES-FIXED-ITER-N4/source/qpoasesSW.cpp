#include "qpoasesSW.h"

void qpoasesSW (elem H[N_QP *N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP], qpOASES::QProblem *example, int nWSR, int *started, tictoc *SOLVERTimer){
    
    USING_NAMESPACE_QPOASES

    // Si el codigo aun no ha realizado una iteracion, se cargan los todos los datos
    // a las estructuras de datos al solver
    if(*started == 0){
      SOLVERTimer->tic();
      // Argumentos: (Matriz hessiana H, vector funcional o lineal h, matriz de restricciones M, vector de restricciones c, numero de iteraciones)
      example->init( H, h, M, NULL, NULL, NULL, c, nWSR, 0);	// Resuelve el primer problema definiendo todas las matrices
      SOLVERTimer->toc();
      *started = 1; 	// Marca que ya se ha hecho una iteracion
    }
    // Si el codigo ya ha hecho una iteracion, al solver solo se cargan los datos 
    // de los vectores h y c
    else if(*started == 1){
      SOLVERTimer->tic();
      // Argumentos: (Vector funcional o lineal h, vector de restricciones c, numero de iteraciones)
      example->hotstart(h, NULL, NULL, NULL, c, nWSR, 0);	// Resuelve el problema actualizando vectores h y c
      SOLVERTimer->toc();
    }

    // Obtener la solucion del problema
    real_t uOpt[N_QP];
    example->getPrimalSolution( uOpt );

    // Se almacena la solucion obtenida 
    for (int i = 0; i<N_QP; i++){
        tk[i] = uOpt[i];
    }
}
