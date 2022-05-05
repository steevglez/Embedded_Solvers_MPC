/**
 * This file is generated by FiOrdOs, a program licensed under GPL
 * by copyright holder Automatic Control Laboratory, ETH Zurich.
 * 
 * If you are interested in using this file commercially,
 * please contact the copyright holder.
 */

#ifndef FIORDOS5X5_CPARAM_SOLVER_H
#define FIORDOS5X5_CPARAM_SOLVER_H

#ifdef USE_REALTYPE_SINGLE
  typedef float realtype;
#else
  typedef double realtype;
#endif

/* <<< struct Params >>> */
#if defined(USE_PARAMSMODE_PTR)  && !(defined(MATLAB_MEX_FILE) || defined(RT))
#  define PARAMS_MACRO(name,len) const realtype *name;
#else /* the default */
#  define PARAMS_MACRO(name,len) realtype name[len];
#endif
typedef struct {
    PARAMS_MACRO(H,25)
    PARAMS_MACRO(g,5)
    PARAMS_MACRO(Ai,150)
    PARAMS_MACRO(bi,30)
} fiordos5x5_Cparam_Params;
#undef PARAMS_MACRO

/* <<< struct Settings >>> */
typedef struct {
    struct {
        int warmstartInner;
    } approach;
    struct {
        realtype init[5];
        int maxit;
    } algoInner;
    struct {
        realtype init[30];
        int maxit;
    } algoOuter;
} fiordos5x5_Cparam_Settings;

/* <<< struct Result >>> */
typedef struct {
    realtype la[30];
    realtype x[5];
    realtype d;
    int iter;
    int exitflag;
} fiordos5x5_Cparam_Result;

/* <<< struct Work >>> */
typedef struct {
    struct {
        realtype inner_g[5];
        realtype inner_c[1];
        realtype traceH0[1];
    } Prob;
    struct {
        realtype init[5];
        realtype glob_LbL;
        realtype glob_UbL;
        realtype tstep;
        int btnum;
        realtype tstepmin;
        int btstopped;
        realtype res_z[5];
        int res_stopcode;
        int res_iter;
    } algoInner;
    struct {
        realtype init[30];
        realtype glob_LbL;
        realtype glob_UbL;
        realtype tstep;
        int btnum;
        realtype tstepmin;
        int btstopped;
        realtype res_z[30];
        int res_stopcode;
        int res_iter;
    } algoOuter;
} fiordos5x5_Cparam_Work;

/* <<< solver INTERFACE >>> */
void fiordos5x5_Cparam_init(fiordos5x5_Cparam_Params *params, fiordos5x5_Cparam_Settings *settings, fiordos5x5_Cparam_Result *result, fiordos5x5_Cparam_Work *work);
void fiordos5x5_Cparam_solve(fiordos5x5_Cparam_Params *params, fiordos5x5_Cparam_Settings *settings, fiordos5x5_Cparam_Result *result, fiordos5x5_Cparam_Work *work);
#endif
