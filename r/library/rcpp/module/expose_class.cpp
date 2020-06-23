#include <Rcpp.h>

using namespace Rcpp;

class Uniform {
public:
    Uniform(double min_, double max_) : min(min_), max(max_) {}

    NumericVector draw(int n) {
        RNGScope scope;
        return runif(n, min, max);
    }

private:
    double min, max;
};

/// create an external pointer to a Uniform object
RcppExport SEXP Uniform__new(SEXP min_, SEXP max_) {
    // convert inputs to appropriate C++ types
    double min = as<double>(min_), max = as<double>(max_);

    // create a pointer to an Uniform object and wrap it as an external pointer
    Rcpp::XPtr<Uniform> ptr(new Uniform(min, max), true);

    // return the external pointer to the R side
    return ptr;
}

/// invoke the draw method
RcppExport SEXP Uniform__draw(SEXP xp, SEXP n_) {
    // grab the object as a XPtr (smart pointer) to Uniform
    Rcpp::XPtr<Uniform> ptr(xp);
    // convert the parameter to int
    int n = as<int>(n_);

    // invoke the function
    NumericVector res = ptr->draw(n);

    // return the result to R
    return res;
}
