#include <RcppCommon.h>

class Foo {
public:
    Foo();

    // this constructor enables implicit Rcpp::as
    Foo(SEXP);

    // this operator enables implicit Rcpp::wrap
    operator SEXP();
}

#include <Rcpp.h>
