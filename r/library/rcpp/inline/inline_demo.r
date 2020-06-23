library(Rcpp)
library(inline)

incltxt <- "
int fibonacci(const int x) {
    if (x == 0) return(0);
    if (x == 1) return(1);
    return fibonacci(x - 1) + fibonacci(x - 2);
}"

fibRcpp <- cxxfunction(signature(xs = "int"), plugin = "Rcpp", incl = incltxt, body = "
                       int x = Rcpp::as<int>(xs);
                       return Rcpp::wrap( fibonacci(x) );
                       ")

## memoization using C++
mincltxt <- "
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
class Fib {
    public:
        Fib(unsigned int n = 1000) {
            memo.resize(n);
            std::fill( memo.begin(), memo.end(), NAN ); // set to NaN
            memo[0] = 0.0;
            memo[1] = 1.0;
        }
    double fibonacci(int x) {
        if (x < 0) // guard against bad input
            return( (double) NAN );
        if (x >= (int) memo.size())
            throw std::range_error(\"x too large for implementation\");
        if (! std::isnan(memo[x]))
            return(memo[x]);

        // if exist, reuse values
        // build precomputed value via recursion
        memo[x] = fibonacci(x-2) + fibonacci(x-1);
        return( memo[x] );
    }
private:
    std::vector< double > memo;
};
"

mfibRcpp <- cxxfunction(signature(xs = "int"), plugin = "Rcpp", includes = mincltxt, 
    body = "
int x = Rcpp::as<int>(xs);
Fib f;
return Rcpp::wrap( f.fibonacci(x-1) );
")

## linear / iterative solution
fibRcppIter <- cxxfunction(signature(xs = "int"), plugin = "Rcpp", body = "
                           int n = Rcpp::as<int>(xs);
                           double first = 0;
                           double second = 1;
                           double third = 0;
                           for (int i=0; i<n; i++) {
                               third = first + second;
                               first = second;
                               second = third;
                           }
                           return Rcpp::wrap(first);
                           ")


fibRcpp(5)
mfibRcpp(5)
fibRcppIter(5)
