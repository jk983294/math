src <- "
    Rcpp::NumericVector xa(a);
    Rcpp::NumericVector xb(b);
    int n_xa = xa.size(), n_xb = xb.size();

    Rcpp::NumericVector xab(n_xa + n_xb - 1);
    for (int i = 0; i < n_xa; i++)
        for (int j = 0; j < n_xb; j++)
            xab[i + j] += xa[i] * xb[j];
    return xab;
"

fun <- cxxfunction(signature(a = "numeric", b = "numeric"), src, plugin = "Rcpp", 
    verbose = TRUE)
fun(1:4, 2:5)
