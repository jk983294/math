src <- "
    Rcpp::NumericVector invec(vx);
    Rcpp::NumericVector outvec(vx);

    for (int i=0; i<invec.size(); i++) {
        outvec[i] = log(invec[i]);
    }
    return outvec;
"

# invec and outvec share the same underlying pointer to R object
shallow_copy <- cxxfunction(signature(vx = "numeric"), src, plugin = "Rcpp")
x <- seq(1, 3, by = 1)
cbind(x, shallow_copy(x))

src <- "
    Rcpp::NumericVector invec(vx);
    Rcpp::NumericVector outvec = Rcpp::clone(vx);

    for (int i=0; i<invec.size(); i++) {
        outvec[i] = log(invec[i]);
    }
    return outvec;
"

# clone method allocates memory for a new object
deep_copy <- cxxfunction(signature(vx = "numeric"), src, plugin = "Rcpp")
x <- seq(1, 3, by = 1)
cbind(x, deep_copy(x))

# sugar
src <- "
    Rcpp::NumericVector invec(vx);
    Rcpp::NumericVector outvec = log(invec);
    return outvec;
"

deep_copy_sugar <- cxxfunction(signature(vx = "numeric"), src, plugin = "Rcpp")
x <- seq(1, 3, by = 1)
cbind(x, deep_copy_sugar(x))
