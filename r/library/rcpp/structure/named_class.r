someVec <- c(mean = 1.23, dim = 42, cnt = 12)
someVec


src <- "
    Rcpp::NumericVector x =
    Rcpp::NumericVector::create(
            Rcpp::Named(\"mean\") = 1.23,
            Rcpp::Named(\"dim\") = 42,
            Rcpp::Named(\"cnt\") = 12);
    return x; "
fun_create <- cxxfunction(signature(), src, plugin = "Rcpp")
fun_create()

# sugar
src <- "
    Rcpp::NumericVector x =
    Rcpp::NumericVector::create(
            _[\"mean\"] = 1.23,
            _[\"dim\"] = 42,
            _[\"cnt\"] = 12);
    return x; "
fun_create_sugar <- cxxfunction(signature(), src, plugin = "Rcpp")
fun_create_sugar()
