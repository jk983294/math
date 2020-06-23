#include <Rcpp.h>

/**
 * demo three way to use R math lib in cpp file
 */
extern "C" SEXP mypnorm(SEXP xx) {
    Rcpp::NumericVector x(xx);
    int n = x.size();
    Rcpp::NumericVector y1(n), y2(n), y3(n);
    for (int i = 0; i < n; i++) {
        // accessing function via remapped R header
        y1[i] = ::Rf_pnorm5(x[i], 0.0, 1.0, 1, 0);

        // or accessing same function via Rcpp’s ’namespace R’
        y2[i] = R::pnorm(x[i], 0.0, 1.0, 1, 0);
    }
    // or using Rcpp sugar which is vectorized
    y3 = Rcpp::pnorm(x);

    return Rcpp::DataFrame::create(Rcpp::Named("R") = y1, Rcpp::Named("Rf_") = y2, Rcpp::Named("sugar") = y3);
}
