#include <iostream>
#include <fstream>
#include "specs.h"
#include <math.h>
#include "mpcSW.h"
#include "tictoc.h"
#include <iomanip>

using namespace std;
//#define DISPLAY

int main(int argc, char *argv[]){

    int OUTPUT = atoi(argv[2]);
    int TIME = atoi(argv[3]);
    int NUM_ITER = stoi(argv[5]);
	
    cout << std::setprecision(15);
    cout << "MPC testbench" << endl;

    elem x[N_SYS];
    elem yref[P_SYS];
    elem xmin[N_SYS];
    elem xmax[N_SYS];
    elem umin[M_SYS];
    elem umax[M_SYS];
    elem A[N_SYS][N_SYS];
    elem B[N_SYS];
    elem Acal[N_SYS*N_QP][N_SYS];
    elem AcalOmgOcal[N_SYS][N_QP];
    elem H[N_QP *N_QP];
    elem Mx[I_QP *N_QP];
    elem L_invLast[(N_SYS+P_SYS)];
    elem expectedResult[P_SYS];
    elem u[P_SYS];
    int errors = 0;

    // Indica problemas con los argumentos del programa
    if (argc!=6){
        cerr << "Must specify .txt\n";
        return 1;
    }

    // Chequear si el archivo de muestras abre correctamente
    ifstream samples(argv[1], ios::in);
    if (!samples.is_open()) {
        cerr << "There was a problem opening the input file: ";
        cerr << argv[1] << endl;
        return 1;
    }

    // Umbral para la diferencia entre los resultados obtenidos y los de referencia(MATLAB)
    elem threshold = 1e-4;

    int nSamples = 0;
    int samplesI = 0;
    int samplesN = 0;

    // Timers para lazo MPC y solver
    tictoc MPCTimer(nSamples);
    tictoc SOLVERTimer(nSamples);
    
    // Obtener numero de dimension del sistema(horizonte de prediccion), restricciones y muestras
    string line;
    getline(samples, line);
    sscanf(line.c_str(), "N:\t%d\tI: %d\tnSamples:\t%d", &samplesN, &samplesI, &nSamples);

    // Chequear si N e I corresponden a lo indicado en N_QP e I_QP en specs.h
    if (N_QP != samplesN){
        cerr << "N_QP from " << argv[1] << " (" << samplesN << ") does not match N_QP from specs.h (";
        cerr << N_QP << ")" << endl;
        return 1;
    }
    if (I_QP != samplesI){
        cerr << "I_QP from " << argv[1] << " (" << samplesI << ")does not match I_QP from specs.h (";
        cerr << I_QP << ")" << endl;
        return 1;
    }

    // Cargar vector xmin
    for (int r=0; r<N_SYS; r++){
    	samples >> xmin[r];
    }
    // Cargar vector xmax
    for (int r=0; r<N_SYS; r++){
    	samples >> xmax[r];
    }
    // Cargar vector umin
    for (int r=0; r<M_SYS; r++){
    	samples >> umin[r];
    }
    // Cargar vector umax
    for (int r=0; r<M_SYS; r++){
    	samples >> umax[r];
    }
    // Cargar matriz A del SISTEMA FISICO (No confundir con matriz de restricciones)
    for (int r=0; r<N_SYS; r++){
    	for (int c=0; c<N_SYS; c++){
    		samples >> A[r][c];
    	}
    }
    // Cargar matriz B del SISTEMA FISICO
    for (int r=0; r<N_SYS; r++){
    	samples >> B[r];
    }
    // Cargar matriz Acal
    for (int r=0; r<N_SYS*N_QP; r++){
    	for (int c=0; c<N_SYS; c++){
    		samples >> Acal[r][c];
    	}
    }
    // Cargar matriz AcalOmgOcal
    for (int r=0; r<N_SYS; r++){
    	for (int c=0; c<N_QP; c++){
    		samples >> AcalOmgOcal[r][c];
    	}
    }
    // Cargar matriz hessiana H
    for (int r=0; r<N_QP; r++){
    	for (int c=0; c<N_QP; c++){
    		samples >> H[r*N_QP +c];
    	}
    }
    // Cargar matriz de restricciones Mx
    for (int r=0; r<I_QP; r++){
        for (int c=0; c<N_QP; c++){
	     samples >> Mx[r*N_QP +c];
		}
    }
    // Cargar vector L_last
	for (int r=0; r<N_SYS+P_SYS; r++){
		samples >> L_invLast[r];
    }
    
    // Iniciacion de variables de estado del sistema
    elem x0[2] = {0, 0};
    elem x_real[2];
    x_real[0] = x0[0];
    x_real[1] = x0[1];
    
    double output_array[3*nSamples];
    
    // Definicion e iniciacion de estructuras del solver
    OSQPWorkspace *work;
    OSQPSettings  *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
    OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));
    
    // Configuracion del solver
    if (settings){
    	osqp_set_default_settings(settings);
    	settings->verbose = 0;
    	settings->max_iter = NUM_ITER;
    }
    
    // Iniciacion de variable que indica si ya inicio una iteracion
    int started = 0;

    for (int sample=0; sample<nSamples; sample++){
    	// Terminar si se acaba el archivo de muestras
        if (samples.peek() == EOF) break;
        // Cargar estados - SOLO DISPONIBLE PARAREFERENCIA, NO USADO EN CODIGO NI LAZO
        for (int i=0; i<N_SYS; i++){
            samples >> x[i];
        }
        // Cargar vector de referencia
        for (int i=0; i<P_SYS; i++){
            samples >> yref[i];
        }
        // Cargar vector de resultados esperados - SOLO REFERENCIA, NO USADO EN EL LAZO
        for (int i=0; i<P_SYS; i++){
            samples >> expectedResult[i];
        }

	// Funcion de lazo MPC
        MPCTimer.tic();
        mpcSW(x_real, yref, xmin, xmax, umin, umax, A, B, Acal, AcalOmgOcal, H, Mx, L_invLast, u,  work, settings, data, &started, output_array, &SOLVERTimer, sample);
        MPCTimer.toc();

        // Calculo de diferencias entre la actuacion obtenida y la de referencia, solo como primera idea de los resultados
        for (int i=0; i<P_SYS; i++){
        	elem error = fabs(expectedResult[i] -u[i]);
            if ((error > threshold) | (error != error)){
                errors++;
        // Mostrar diferencias
#ifdef DISPLAY
                cout << "sample number : " << sample << " " << i << endl;
                cout << expectedResult[i] << "\t";
                cout << u[i] << "\t";
                cout << error;
                cout << "\n";
#endif

            }
        }

    }
    
    // Limpieza de estructuras del solver
    osqp_cleanup(work);
    if (data) {
        if (data->A) c_free(data->A);
        if (data->P) c_free(data->P);
        c_free(data);
    }
    if (settings)  c_free(settings);
    
    char dim[5];
    
    if(OUTPUT == 1){
    	// CREACION Y ESCRITURA DE ARCHIVO .CSV
	// Nombre de archivo csv
	sprintf(dim, "%d", N_QP);
        char csv_name_out[35];
        char format_csv_name_out[] = "csv_osqpMPC%sx%sOUT_iter%d.csv";
        sprintf(csv_name_out, format_csv_name_out, dim, dim, NUM_ITER);
	// Creacion del archivo y escritura al final del archivo
	FILE *fp;
	fp = fopen(csv_name_out, "a");
	// Imprimir
	for(int i = 0 ; i < nSamples ; i++)
	    fprintf(fp, "%.7f,%.7f,%.7f\n", output_array[3*i], output_array[3*i+1], output_array[3*i+2] );
	// Cierre del archivo y termino
	fclose(fp);
    }
    
    if(TIME == 1){
	// CREACION Y ESCRITURA DE ARCHIVO .CSV
	// Nombre de archivo csv
	sprintf(dim, "%d", N_QP);
    	char csv_name_time[40];
    	char format_csv_name_time[] = "csv_osqpMPC%sx%sTIME_iter%d_%s.csv";
    	sprintf(csv_name_time, format_csv_name_time, dim, dim, NUM_ITER, argv[4]);
	// Creacion del archivo y escritura al final del archivo
	FILE *fp2;
	fp2 = fopen(csv_name_time, "a");
	// Escribir
	for(int i = 0 ; i < nSamples ; i++)
	    fprintf(fp2, "%.3f,%.3f\n", MPCTimer.tocs[i]/1e3, SOLVERTimer.tocs[i]/1e3);
	// Cierre del archivo y termino
	fclose(fp2);
    }

    samples.close();

    // Imprimir informacion importante sobre muestras y tiempos
    cout << "Finished processing " << nSamples << " samples" << endl;
    cout << "Number of differences between expected and calculated:\t" << errors << endl;
    cout << "Threshold: " << threshold << endl;
    std::cout << std::endl;
    double sw_cycles_mpc = MPCTimer.mean();
    std::cout << "MPC TIMER" << std::endl;
    std::cout << "Mean time running application in software: \t\t" << sw_cycles_mpc / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Standard deviation running application in software: \t" << MPCTimer.stdev() / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Max running application in software: \t\t\t" << MPCTimer.max() / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Min running application in software: \t\t\t" << MPCTimer.min() / 1e3 << "\xC2\xB5s" << std::endl;
    double sw_cycles_solver = SOLVERTimer.mean();
    std::cout << "SOLVER TIMER" << std::endl;
    std::cout << "Mean time running application in software: \t\t" << sw_cycles_solver / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Standard deviation running application in software: \t" << SOLVERTimer.stdev() / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Max running application in software: \t\t\t" << SOLVERTimer.max() / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Min running application in software: \t\t\t" << SOLVERTimer.min() / 1e3 << "\xC2\xB5s" << std::endl;

    return errors;
}
