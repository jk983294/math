src <- "
    Rcpp::NumericMatrix mat = Rcpp::clone<Rcpp::NumericMatrix>(mx);
    std::transform(mat.begin(), mat.end(), mat.begin(), ::sqrt);
    return mat;
"

# Internally, matrices are implemented as vectors with an associated dimension
# attribute
fun <- cxxfunction(signature(mx = "numeric"), src, plugin = "Rcpp")
orig <- matrix(1:9, 3, 3)
fun(orig)
