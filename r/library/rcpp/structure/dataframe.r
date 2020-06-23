# Internally, data frames are represented as lists different columns will always
# be recycled to have common length.

src <- "
    Rcpp::IntegerVector v = Rcpp::IntegerVector::create(7,8,9);
    std::vector<std::string> s(3);
    s[0] = \"x\";
    s[1] = \"y\";
    s[2] = \"z\";
    return Rcpp::DataFrame::create(Rcpp::Named(\"a\")=v, Rcpp::Named(\"b\")=s);
"

fun_create <- cxxfunction(signature(), src, plugin = "Rcpp")
fun_create()
