// Compilar: g++ osqp5x5_1.cpp -IOSQP_FOLDER/include OSQP_FOLDER/lib/libosqp.a -o osqp5x5_1
// Ejecutar: ./osqp5x5_1

#include "osqpSW.h"

void H_format_OSQP(double M[], c_float values[], c_int index[], c_int column_pointer[], c_int * P_nnz);
void M_format_OSQP(double M[], c_float values[], c_int index[], c_int column_pointer[], c_int * A_nnz, c_float l[]);

void osqpSW(elem H[N_QP*N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP], OSQPWorkspace * & work, OSQPSettings *settings, OSQPData *data,  int *started, tictoc *SOLVERTimer){

    
    // Si el codigo aun no ha realizado una iteracion, se cargan los todos los datos
    // a las estructuras de datos al solver
    if(*started == 0){
    
        c_int n = N_QP;	// Dimension de matriz M es nxn, segun el horizonte de prediccion
    	c_int m = I_QP;	// Dimension de matriz A es mxn, segun el horizonte de prediccion y numero de restricciones
    	
    	// Matriz H - Hessiana en formato Compresed Sparse Column(CSC) (Esta matriz se presenta con la letra P en OSQP)
    	c_float P_x[(N_QP*(N_QP+1))/2];		// Elementos no nulos en M
    	c_int   P_i[(N_QP*(N_QP+1))/2];		// Indice de fila de elementos no nulos en M
    	c_int   P_p[(N_QP+1)];				// Puntero de columnas a elementos no nulos en M
    	c_int   P_nnz;					// Numero de elementos no nulos en M
    	H_format_OSQP(H, P_x, P_i, P_p, &P_nnz);	// Funcion que convierte M en formato CSC
    
    	// Matriz M(Mx) y vector l - Matriz de restriciones M en formato Compresed Sparse Column(CSC) (Esta matriz se presenta con la letra A en OSQP)
    	c_float A_x[N_QP*I_QP];			// Elementos no nulos en M
    	c_int A_i[N_QP*I_QP];				// Indice de fila de elementos no nulos en M
    	c_int A_p[N_QP+1];				// Puntero de columnas a elementos no nulos en M
    	c_int A_nnz;					// Numero de elementos no nulos en M
    	c_float l[I_QP];				// Vector l
    	M_format_OSQP(M, A_x, A_i, A_p, &A_nnz, l);	// Funcion que convierte A en formato CSC y ajusta l

    	// Exitflag
    	c_int exitflag = 0;

    	// Cargar vectores a la estructura de datos data de OSQP
    	if (data) {
            data->n = n; // dimension n	
            data->m = m; // dimension m
            data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p); // Matriz H o P - Hessiana
            data->q = h;						     // Vector h o q - Funcional o lineal (Este vector se presenta con la letra q en OSQP)
            data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p); // Matriz M o A - Matriz de restricciones
            data->l = l;						     // Vector l - No utilizado, pero si definido a -inf			
            data->u = c;						     // Vector c - Vector de restricciones (Este vector se presenta con la letra u en OSQP)
        }

        // Configuracion del workspace de OSQP
        exitflag = osqp_setup(&work, data, settings);
        
        *started = 1;		// Marca que ya se ha hecho una iteracion        
    }
    // Si el codigo ya ha hecho una iteracion, al solver solo se cargan los datos 
    // de las matrices h y c
    else if(*started == 1){
    	osqp_update_lin_cost(work, h); 	// Actualizacion de vector h o q
	osqp_update_upper_bound(work, c); 	// Actualizacion de vector c o u
    }

    // Resolver el problema
    SOLVERTimer->tic();
    osqp_solve(work);
    SOLVERTimer->toc();

    // Se almacena la solucion obtenida 
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
    // Vector l no es usado, pero debe definirse. Por lo que se ajusta al menor valor posible proporcionado por OSQP
    for(i = 0 ; i < I_QP ; i++){
      l[i] = -OSQP_INFTY;
    }
}

