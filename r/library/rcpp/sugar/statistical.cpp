#include <Rcpp.h>

x1 = dnorm(y1, 0, 1);  // density of y1 at m=0, sd=1
x2 = pnorm(y2, 0, 1);  // distribution function of y2
x3 = qnorm(y3, 0, 1);  // quantiles of y3
x4 = rnorm(n, 0, 1);   // ’n’ RNG draws of N(0, 1)
