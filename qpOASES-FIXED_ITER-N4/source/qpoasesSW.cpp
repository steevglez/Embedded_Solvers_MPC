// Instrucciones de programa anterior:
// Compilar: g++ qpoases5x5_1.cpp -IqpOASES_FOLDER/include qpOASES_FOLDER/build/libs/libqpOASES.a -o qpoases5x5_1
// Ejecutar: ./qpoases5x5_1

#include "qpoasesSW.h"

void qpoasesSW (elem H[N_QP *N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP], qpOASES::QProblem *example, int nWSR, int *started, tictoc *SOLVERTimer){
    USING_NAMESPACE_QPOASES

    if(*started == 0){
      SOLVERTimer->tic();
      example->init( H, h, M, NULL, NULL, NULL, c, nWSR,0 );
      SOLVERTimer->toc();
      *started = 1;
    }else if(*started == 1){
      SOLVERTimer->tic();
      example->hotstart(h, NULL, NULL, NULL, c, nWSR, 0);
      SOLVERTimer->toc();
    }

    /* Get and print solution of first QP. */
    real_t uOpt[N_QP];
    example->getPrimalSolution( uOpt );

    for (int i = 0; i<N_QP; i++){
        tk[i] = uOpt[i];
    }
    // uOpt[0], uOpt[1], uOpt[2], uOpt[3], uOpt[4], example.getObjVal()
}
