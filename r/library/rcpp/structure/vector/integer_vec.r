library(inline)

src <- "
    Rcpp::IntegerVector epn(4);
    epn[0] = 6;
    epn[1] = 14;
    epn[2] = 496;
    epn[3] = 8182;
    return epn;
"

fun_output_vector <- cxxfunction(signature(), src, plugin = "Rcpp")
fun_output_vector()


src <- "
    Rcpp::IntegerVector vec(vx);
    int prod = 1;
    for (int i=0; i<vec.size(); i++) {
        prod *= vec[i];
    }
    return Rcpp::wrap(prod);
"
fun_input_vector <- cxxfunction(signature(vx = "integer"), src, plugin = "Rcpp")
fun_input_vector(1:10)


# wrong inputs
fun_input_vector(seq(1, 1.9, by = 0.1))
fun_input_vector(LETTERS[1:10])
