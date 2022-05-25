#ifndef MCMCHH
#define MCMCHH

#include "mcmc.h"
#include "mcmc.hpp"

enum {
    nburn = 500,
    ndraw = 10000
};

void genmcmc(double lambda[], arma::mat &draws_out); // lambda size if NBins

#endif
