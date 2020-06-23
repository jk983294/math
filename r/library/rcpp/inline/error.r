library(inline)

src <- "
    int dx = Rcpp::as<int>(x);
    if( dx > 10 )
        throw std::range_error(\"too big\");
    return Rcpp::wrap( dx * dx);
"
fun <- cxxfunction(signature(x = "integer"), body = src, plugin = "Rcpp")
fun(3)
fun(13)
