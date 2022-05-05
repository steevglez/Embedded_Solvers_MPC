#ifndef SPECS_H
#define SPECS_H

#ifndef N_QP
#define N_QP 4 	// Dimension del sistema segun el valor del horizonte de prediccion N - Modificable
#endif
#ifndef I_QP
#define I_QP (6*N_QP)	// Numero de restricciones - No modificar a menos que el modelo fisico cambie
#endif

#define N_SYS 2
#define M_SYS 1
#define P_SYS 1

// Definine en tipo de dato, double o single. Especificar solo uno
#define ELEM_DOUBLE
//#define ELEM_FLOAT

#if defined ELEM_DOUBLE
typedef double elem;
#elif defined ELEM_FLOAT
typedef float elem;
#else
#error Something is wrong with how elem is defined
#endif

#endif
