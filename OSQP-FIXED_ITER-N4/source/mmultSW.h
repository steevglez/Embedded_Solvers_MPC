#ifndef MMULTSW_H
#define MMULTSW_H

#include "specs.h"
#include <cblas.h>

void mmultSW (elem *A, CBLAS_TRANSPOSE tansposeA, elem *B, CBLAS_TRANSPOSE tansposeB, elem *C, int m, int n, int p);

#endif
