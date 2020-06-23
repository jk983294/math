#include <Rcpp.h>

int fibonacci(const int x) {
    if (x == 0) return (0);
    if (x == 1) return (1);
    return (fibonacci(x - 1)) + fibonacci(x - 2);
}

/**
 * SEXP = pointers to S expression objects
 * Essentially everything inside R is represented as such a SEXP object
 */
extern "C" SEXP fibWrapper(SEXP xs) {
    int x = Rcpp::as<int>(xs);
    int fib = fibonacci(x);
    return (Rcpp::wrap(fib));
}

extern "C" SEXP fun(SEXP x) {
    try {
        int dx = Rcpp::as<int>(x);
        if (dx > 10) throw std::range_error("too big");
        return Rcpp::wrap(dx * dx);
    } catch (std::exception& __ex__) {
        forward_exception_to_r(__ex__);
    } catch (...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;  // not reached
}
