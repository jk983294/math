library(inline)

src <- "
    Rcpp::NumericVector a(4);
    a[0] = 6.0;
    a[1] = 14.0;
    a[2] = 496.0;
    a[3] = 8182.0;
    return a;
"

fun_output_vector <- cxxfunction(signature(), src, plugin = "Rcpp")
fun_output_vector()


src <- "
    Rcpp::NumericVector vec(vx);
    double p = Rcpp::as<double>(dd);
    double sum = 0;
    for (int i=0; i<vec.size(); i++) {
        sum += pow(vec[i], p);
    }
    return Rcpp::wrap(sum);
"
fun_input_vector <- cxxfunction(signature(vx = "numeric", dd = "numeric"), src, plugin = "Rcpp")
fun_input_vector(1:4, 2)
fun_input_vector(1:4, 2.2)


# wrong inputs
fun_input_vector(seq(1, 1.9, by = 0.1))
fun_input_vector(LETTERS[1:10])
