list_data <- list("Red", c(1.5, 2.5), il = c(1, 2), bv = TRUE, 51.23)
list_data

src <- "
    Rcpp::List control(l);
    std::string str = Rcpp::as<std::string>(control[0]);
    Rcpp::NumericVector dl = Rcpp::as<Rcpp::NumericVector>(control[1]);
    Rcpp::IntegerVector il = Rcpp::as<Rcpp::IntegerVector>(control[\"il\"]);
    bool bv = Rcpp::as<bool>(control[\"bv\"]);
    double dv = Rcpp::as<double>(control[4]);

    return Rcpp::List::create(Rcpp::Named(\"str\") = str,
                            Rcpp::Named(\"dl\") = dl,
                            Rcpp::Named(\"il\") = il,
                            Rcpp::Named(\"bv\") = bv,
                            Rcpp::Named(\"dv\") = dv);
"
fun <- cxxfunction(signature(l = "list"), src, plugin = "Rcpp")
fun(list_data)
