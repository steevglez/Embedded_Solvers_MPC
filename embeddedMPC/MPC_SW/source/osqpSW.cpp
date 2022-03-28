// Compilar: g++ osqp5x5_1.cpp -IOSQP_FOLDER/include OSQP_FOLDER/lib/libosqp.a -o osqp5x5_1
// Ejecutar: ./osqp5x5_1

#include "osqpSW.h"

void osqpSW(elem H[N_QP *N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP]){

    // Load problem data

    // H
    // CSC_vectors(H)

    c_float P_x[15] = {0.2008398256540883, 0.0008380730958730651, 0.20083637427514697, 0.0008362739504127661, 0.0008346281186590991, 0.20083293568778623, 0.0008344298698261268, 0.0008328355561809345, 0.0008311961020701798, 0.2008295102294131, 0.0008325424712014894, 0.0008309982449631417, 0.0008294102964048211, 0.0008277773875982911, 0.20082609824556535};
    c_int   P_nnz  = 15;
    c_int   P_i[15] = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4};
    c_int   P_p[6] = {0, 1, 3, 6, 10, 15};

    c_float q[5]   = {-0.6442398069487796, -0.6416124969399188, -0.6389900601491728, -0.6363726345539936, -0.6337603620384462};
    // c_float q[5] = h


    // M
    // CSC_vectors(M)

    c_float A_x[70] = {0.038539485272602474, 1.935941000707434e-05, 0.03747834528928226, 5.73658564010647e-05, 0.03644642256343802, 9.432583941840871e-05, 0.03544291263180586, 0.00013026817221325687, 0.034467033181083906, 0.0001652208746030117, -0.038539485272602474, -1.935941000707434e-05, -0.03747834528928226, -5.73658564010647e-05, -0.03644642256343802, -9.432583941840871e-05, -0.03544291263180586, -0.00013026817221325687, -0.034467033181083906, -0.0001652208746030117, 1.0, -1.0, 0.038539485272602474, 1.935941000707434e-05, 0.03747834528928226, 5.73658564010647e-05, 0.03644642256343802, 9.432583941840871e-05, 0.03544291263180586, 0.00013026817221325687, -0.038539485272602474, -1.935941000707434e-05, -0.03747834528928226, -5.73658564010647e-05, -0.03644642256343802, -9.432583941840871e-05, -0.03544291263180586, -0.00013026817221325687, 1.0, -1.0, 0.038539485272602474, 1.935941000707434e-05, 0.03747834528928226, 5.73658564010647e-05, 0.03644642256343802, 9.432583941840871e-05, -0.038539485272602474, -1.935941000707434e-05, -0.03747834528928226, -5.73658564010647e-05, -0.03644642256343802, -9.432583941840871e-05, 1.0, -1.0, 0.038539485272602474, 1.935941000707434e-05, 0.03747834528928226, 5.73658564010647e-05, -0.038539485272602474, -1.935941000707434e-05, -0.03747834528928226, -5.73658564010647e-05, 1.0, -1.0, 0.038539485272602474, 1.935941000707434e-05, -0.038539485272602474, -1.935941000707434e-05, 1.0, -1.0};
    c_int   A_nnz  = 70;
    c_int   A_i[70] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 17, 18, 19, 21, 26, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 22, 27, 6, 7, 8, 9, 16, 17, 18, 19, 23, 28, 8, 9, 18, 19, 24, 29};
    c_int   A_p[6] = {0, 21, 39, 53, 63, 70};

    c_float l[30]   = {-OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY};
    c_float u[30]   = {5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    // c_float u[24] = c

    c_int n = 5;
    c_int m = 30;

    // Exitflag
    c_int exitflag = 0;

    // Workspace structures
    OSQPWorkspace *work;
    OSQPSettings  *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
    OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));

    // Populate data
    if (data) {
        data->n = n;
        data->m = m;
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
        data->q = q;
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
        data->l = l;
        data->u = u;
    }

    // Define solver settings as default
    if (settings) osqp_set_default_settings(settings);

    // Setup workspace
    exitflag = osqp_setup(&work, data, settings);

    // Solve Problem
    osqp_solve(work);

    for(int i = 0 ; i < col ; i++){
        tk[i] = (work->solution->x+i);
    }

    // Clean workspace
    osqp_cleanup(work);
    if (data) {
        if (data->A) c_free(data->A);
        if (data->P) c_free(data->P);
        c_free(data);
    }
    if (settings)  c_free(settings);

}

