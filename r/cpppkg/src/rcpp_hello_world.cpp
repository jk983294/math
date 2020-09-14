#include <Rcpp.h>
#include <math_dummy.h>

using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_hello_world() {
    CharacterVector x = CharacterVector::create("foo", "bar");
    NumericVector y = NumericVector::create(0.0, 1.0);
    List z = List::create(x, y);

    return z;
}

//' Multiply a number by two
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
int MyTimesTwo(int x) { return ornate::dummy_prod(x, 2); }
