#include "pysnow.h"
#include "pySNvC.h"
#include "fisher.h"
#include "mcmc.h"
#include "mcmc.hpp" // mcmc package
#include <cmath>
#include <iostream>

#define expr expr_junoes

using namespace std;


static double default_x[] = {2.5, 2.5, 2.5, 9.5, 12., 15.6, 5., 5., 5.};
int main() {

    double lambda[NBins] = {0.};
    double dlambda[ndim][NBins] = {0.};
    rateGen(expr, default_x, 10., lambda, NBins);
//    pygetSpecArgo(default_x, lambda);
//    pygetSpecResNova(default_x, lambda);
    diff_5p(default_x, 1e-3, dlambda);
    /*
    for(size_t k = 0; k != ndim; ++k) {
        for(size_t l = 0; l != NBins; ++l) {
            printf("%lf ", dlambda[k][l]);
        }
        printf("\n");
    }
    */

    char labels[ndim][10] = {"ae", "aebar", "ax", 
                             "<Ee>", "<Eebar>", "<Ex>",
                             "Ee", "Eebar", "Ex"};
    arma::mat sampler;
    genmcmc(lambda, sampler);
    double N[NBins];
    double tmp, sum = 0.;
    printf("fisher[ndim]:\n");
    for(size_t i = 0; i != ndim; ++i) {
        printf("%10s", labels[i]); // vertical label
        for(size_t j = 0; j <= i; ++j) {
            for(size_t k = 0; k != ndraw; ++k) {
                for(size_t bin = 0; bin != NBins; ++bin) {
                    N[bin] = sampler(k, bin);
                }
                if(i == j) {
                    tmp = pLLHi(N, lambda, dlambda[i]);
                    sum += tmp * tmp;
                } else {
                    sum += pLLHi(N, lambda, dlambda[i]) * pLLHi(N, lambda, dlambda[j]);
                }
            }
            printf("% 10.2E", sum/ndraw);
            sum = 0.;
        }
        printf("\n");
        fflush(stdout);
    }
    printf("%10s", " ");
    for(size_t i = 0; i != ndim; ++i) printf("%10s", labels[i]); printf("\n");
    fflush(stdout);

    return 0;
}

double pLLHi(double N[NBins], double lambda[NBins], double dlambda[NBins]) {
    
    double res = 0.;
    for(size_t i = 0; i != NBins; ++i) {
        res += dlambda[i] * (N[i]/lambda[i] - 1);
    }

    return res;
}

void diff_5p(double x0[], double h, double res[ndim][NBins]) {
    
    double x[ndim];
    for(size_t k = 0; k != ndim; ++k) {
        x[k] = x0[k];
    }
    for(size_t i = 0; i != ndim; ++i) {
        for(size_t j = 0; j != NBins; ++j) {
            res[i][j] = 0.;
        }
    }
    
    double oneover12h = 1./(12*h);
    double lamb[NBins] = {0.};
    double delta[4] = {-2*h, -h, h, 2*h};

    static double coeff[4] = {1, -8, 8, -1};
    for(size_t k = 0; k != ndim; ++k) {
        for(size_t hindex = 0; hindex != 4; ++hindex) {
            x[k] = x0[k] + delta[hindex];
            rateGen(expr, x, 10., lamb, NBins);
//            pygetSpecArgo(x, lamb);
//            pygetSpecResNova(x, lamb);
            for(size_t bin = 0; bin != NBins; ++bin) {
                res[k][bin] += coeff[hindex] * lamb[bin] * oneover12h;
            }
        }
        x[k] = x0[k];
    }
}





