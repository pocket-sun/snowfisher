#include "mcmc.hpp"
#include "fisher.h"
#include "mcmc.h"

double ll_dens(const arma::vec& vals_inp, void* params) {

    double *lambdas = (double*) params;
    double ret = 0.;
    double nr = 1.;

    for(unsigned int i = 0; i != NBins; ++i) {
        nr = std::round(vals_inp(i));
        ret += -lambdas[i] 
         + nr * std::log(lambdas[i])
         - arma::accu(arma::log(arma::linspace(1,nr,nr)));
    }

    return ret;
}

void genmcmc(double lambda[], arma::mat &draws_out) { // lambda size if NBins

    arma::vec initial_val(NBins);
    for(size_t k = 0; k != NBins; ++k) {
        initial_val(k) = lambda[k];
    }

    mcmc::algo_settings_t settings;
    settings.rwmh_n_burnin = nburn;
    settings.rwmh_n_draws = ndraw;
    settings.vals_bound = true;
    arma::vec lb(NBins); 
    arma::vec ub(NBins); 
    for(size_t k = 0; k != NBins; ++k) {
        lb(k) = 1.;
        ub(k) = lambda[k]*100;
    }
    settings.lower_bounds = lb;
    settings.upper_bounds = ub;

    //arma::mat draws_out;
    mcmc::rwmh(initial_val,draws_out,ll_dens,&lambda,settings);
    draws_out = arma::round(draws_out);

//    arma::cout << draws_out.rows(0, 15);

}

