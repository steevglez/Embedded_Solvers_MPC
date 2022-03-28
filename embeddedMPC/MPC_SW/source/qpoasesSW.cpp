// Instrucciones de programa anterior:
// Compilar: g++ qpoases5x5_1.cpp -IqpOASES_FOLDER/include qpOASES_FOLDER/build/libs/libqpOASES.a -o qpoases5x5_1
// Ejecutar: ./qpoases5x5_1

#include "qpoasesSW.h"

void qpoasesSW (elem H[N_QP *N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP]){
    USING_NAMESPACE_QPOASES

    /* Setting up QProblemB object. */
	QProblemB example( 4 );
	Options options;
	//options.enableFlippingBounds = BT_FALSE;
	options.initialStatusBounds = ST_INACTIVE;
	options.numRefinementSteps = 1;
	options.enableCholeskyRefactorisation = 1;
	example.setOptions( options );

    /* Solve first QP. */
	int_t nWSR = 10;
	example.init( H, h, M, c, nWSR,0 );

    /* Get and print solution of first QP. */
    real_t xOpt[4];
	example.getPrimalSolution( xOpt );

    for (int i = 0; i<N_QP; i++){
        tk[i] = xOpt[i];
    }
    // xOpt[0], xOpt[1], xOpt[2], xOpt[3], xOpt[4], example.getObjVal()
}