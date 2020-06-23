#include <Rcpp.h>

// Rcpp sugar defines the usual binary arithmetic operators: +, -, * , /
// two numeric vectors of the same size
NumericVector x;
NumericVector y;

// expressions involving two vectors
NumericVector res = x + y;
NumericVector res = x - y;
NumericVector res = x * y;  // NB element-wise multiplication
NumericVector res = x / y;

// one vector, one single value
NumericVector res = x + 2.0;
NumericVector res = 2.0 - x;
NumericVector res = y * 2.0;
NumericVector res = 2.0 / y;

// two expressions
NumericVector res = x * y + y / 2.0;
NumericVector res = x * (y - 2.0);
NumericVector res = x / (y * y);

// expressions involving two vectors
LogicalVector res = x < y;
LogicalVector res = x > y;
LogicalVector res = x <= y;
LogicalVector res = x >= y;
LogicalVector res = x == y;
LogicalVector res = x != y;

// one vector, one single value
LogicalVector res = x < 2;
LogicalVector res = 2 > x;
LogicalVector res = y <= 2;
LogicalVector res = 2 != y;

// two expressions
LogicalVector res = (x + y) < (x * x);
LogicalVector res = (x + y) >= (x * x);
LogicalVector res = (x + y) == (x * x);

// negate x
NumericVector res = -x;

// use it as part of a numerical expression
NumericVector res = -x * (x + 2.0);

// negate the logical expression "y < x"
LogicalVector res = !(y < x);
