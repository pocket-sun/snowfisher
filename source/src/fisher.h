#ifndef LOGLIKELIHOOD_H
#define LOGLIKELIHOOD_H

#include "fisher.h"
#include "pysnow.h"

enum {
   NBins = NBinsDUNE,
   ndim = 9 // 9 flux parameters
};

void diff_5p(double x0[], double h, double res[ndim][NBins]);
double pLLHi(double N[NBins], double lambda[NBins], double dlambda[NBins]);

#endif
