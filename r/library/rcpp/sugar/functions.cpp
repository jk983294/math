#include <Rcpp.h>

/**
 * Functions Producing a Single Logical Result
 */
IntegerVector x = seq_len(1000);
IntegerVector y = seq_len(1000);

bool res = is_true(all(x < y));
bool res = is_false(any(x < y));
bool res = is_na(any(x < y));

/**
 * Functions Producing Sugar Expressions
 */
IntegerVector x = IntegerVector::create(0, 1, NA_INTEGER, 3);
is_na(x);      // Each element of result evaluates to TRUE if the corresponding input is a missing value
sign(x);       // element whose values are one of 1, 0, −1, or NA, depending on the sign of the input
seq_along(x);  // values go from 1 to the size of the input
IntegerVector x = seq_len(10);  // ith element expands to i
pmin(x, x* x);                  // ith element = the lowest value between the ith element of two expressions
pmax(x* x, 2);                  // ith element = the highest value between the ith element of two expressions
ifelse(x > y, x, 2);            // condition select element
diff(x);                        // res[i] = x[i+1] - x[i]
duplicated(x);                  // res[i] = (x[i] == x[i - 1])
table(x);  // returns a named vector with counts of the occurrences of each element in the input vector

setdiff(x, y);    // the values of x which are not contained in y vector
union_(x, y);     // union of the two vectors
intersect(x, y);  // the intersection of the two vectors
unique(x);        // the subset of unique values among its input vector
sort_unique(x);   // sort(unique( x )) ;

template <typename T>
T square(const T& x) {
    return x * x;
}
sapply(seq_len(10), square<int>);  // applies function to each element to create a new result
lapply(seq_len(10), square<int>);  // similar to sapply except that the result is always a list

template <typename T>
struct sumOfSquares : std::unary_function<T, T> {
    T operator()(const T& x, const T& y) { return x * x + y * y; }
} NumericVector res =
    mapply(seq_len(10), seq_len(10), sumOfSquares<double>());  // mapply permits multiple vectors as input.

int a, b;
clamp(a, x, b);  // values of the vector x between [a, b]
