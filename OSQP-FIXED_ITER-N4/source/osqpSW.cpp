// Compilar: g++ osqp5x5_1.cpp -IOSQP_FOLDER/include OSQP_FOLDER/lib/libosqp.a -o osqp5x5_1
// Ejecutar: ./osqp5x5_1

#include "osqpSW.h"

void H_format_OSQP(double M[], c_float values[], c_int index[], c_int column_pointer[], c_int * P_nnz);
void M_format_OSQP(double M[], c_float values[], c_int index[], c_int column_pointer[], c_int * A_nnz, c_float l[]);

void osqpSW(elem H[N_QP*N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP], OSQPWorkspace * & work, OSQPSettings *settings, OSQPData *data,  int *started, tictoc *SOLVERTimer){

    // Load problem data
    if(*started == 0){
    
        c_int n = N_QP;
    	c_int m = I_QP;
    	// H
    	c_float P_x[(N_QP*(N_QP+1))/2];
    	c_int   P_i[(N_QP*(N_QP+1))/2];
    	c_int   P_p[(N_QP+1)];
    	c_int   P_nnz;
    	H_format_OSQP(H, P_x, P_i, P_p, &P_nnz);
    
    	// M y l
    	c_float A_x[N_QP*I_QP];
    	c_int A_i[N_QP*I_QP];
    	c_int A_p[N_QP+1];
    	c_int A_nnz;
    	c_float l[I_QP];
    	M_format_OSQP(M, A_x, A_i, A_p, &A_nnz, l);

    	// Exitflag
    	c_int exitflag = 0;

    	// Populate data
    	if (data) {
            data->n = n;
            data->m = m;
            data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
            data->q = h;
            data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
            data->l = l;
            data->u = c;
        }

        // Setup workspace
        exitflag = osqp_setup(&work, data, settings);
        
        *started = 1;
        
    }else if(*started == 1){ 
    	osqp_update_lin_cost(work, h);
	osqp_update_upper_bound(work, c);
    }

    // Solve Problem
    SOLVERTimer->tic();
    osqp_solve(work);
    SOLVERTimer->toc();

    for(int i = 0 ; i < N_QP ; i++){
        tk[i] = *(work->solution->x+i);
    }

}

void H_format_OSQP(double M[], c_float values[], c_int index[], c_int column_pointer[], c_int * P_nnz){
    int colptr_check[N_QP];
    int nnz_elements = 0;
    int i,j,x,y;
    for(i = 0 ;  i < N_QP ; i++)
        colptr_check[i] = 0;
    i = 0;
    j = 0;
    for(x = 0 ;  x < N_QP ; x++){
        for(y = 0 ; y < x+1 ; y++){
            if(M[x+N_QP*y] != 0){
                values[i] = M[x+N_QP*y];
                index[i] = y;
                nnz_elements = (nnz_elements + 1);
                i++;
                if(colptr_check[x] == 0){
                    column_pointer[j] = (nnz_elements-1);
                    colptr_check[x] = 1;
                    j++;
                }
            }
        }
    }
    column_pointer[j] = nnz_elements;
    *P_nnz = nnz_elements;
    //values[i] = '\0';
    //index[i] = '\0';
    //column_pointer[j+1] = '\0';
}

void M_format_OSQP(double M[], c_float values[], c_int index[], c_int column_pointer[], c_int * A_nnz, c_float l[]){
    int colptr_check[N_QP];
    int nnz_elements = 0;
    int i,j,x,y;
    for(i = 0 ;  i < N_QP ; i++)
        colptr_check[i] = 0;
    i = 0;
    j = 0;
    for(x = 0 ;  x < N_QP ; x++){
        for(y = 0 ; y < I_QP ; y++){
            if(M[x+N_QP*y] != 0){
                values[i] = M[x+N_QP*y];
                index[i] = y;
                nnz_elements = (nnz_elements + 1);
                i++;
                if(colptr_check[x] == 0){
                    column_pointer[j] = (nnz_elements-1);
                    colptr_check[x] = 1;
                    j++;
                }
            }
        }
    }
    column_pointer[j] = nnz_elements;
    *A_nnz = nnz_elements;
    for(i = 0 ; i < I_QP ; i++){
      l[i] = -OSQP_INFTY;
    }
    //values[i] = '\0';
    //index[i] = '\0';
    //column_pointer[j+1] = '\0';
}

