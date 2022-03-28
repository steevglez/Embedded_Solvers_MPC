#ifndef SPECS_H
#define SPECS_H

// size of system (?)

#ifndef N_QP
#define N_QP 4
#endif
#ifndef I_QP
#define I_QP (6*N_QP)
#endif

#define N_SYS 2
#define M_SYS 1
#define P_SYS 1

// MS_QP is needed for ls_chol()
// the next line is the same as: MS_QP = (int) (ceil((N/2))*ceil((N/2)))
#define MS_QP  ((N_QP/2)+((N_QP%2) != 0)) * ((N_QP/2)+((N_QP%2) != 0))


// type of data (only one at a time)
//#define ELEM_DOUBLE
#define ELEM_FLOAT

#if defined ELEM_DOUBLE
typedef double elem;
#elif defined ELEM_FLOAT
typedef float elem;
#else
#error Something is wrong with how elem is defined
#endif

// linear solver to use (only one at a time)
//#define MINRES
//#define CGWIP
#define CHOL

#endif
