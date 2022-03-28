
#include "mmultSW.h"


void mmultSW(elem *A, CBLAS_TRANSPOSE tansposeA, elem *B, CBLAS_TRANSPOSE tansposeB, elem *C, int m, int n, int p){
    int lda = 0;
    int ldb = 0;

    if (tansposeA == CblasTrans){
        lda = m;
    }
    else if (tansposeA == CblasNoTrans){
        lda = n;
    }
    if (tansposeB == CblasTrans){
        ldb = n;
    }
    else if (tansposeB == CblasNoTrans){
        ldb = p;
    }

#ifdef ELEM_DOUBLE
    cblas_dgemm(CblasRowMajor, tansposeA, tansposeB, m, p, n, 1, A, lda, B, ldb, 0, C, p);
#elif defined ELEM_FLOAT
    cblas_sgemm(CblasRowMajor, tansposeA, tansposeB, m, p, n, 1, A, lda, B, ldb, 0, C, p);
#else
#error Something is wrong with how elem is defined
#endif
}

