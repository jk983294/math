fun <- cxxfunction(signature(), plugin = "Rcpp", body = "
    Rcpp::LogicalVector v(6);
    v[0] = v[1] = false;
    v[1] = true;
    v[3] = R_NaN;
    v[4] = R_PosInf;
    v[5] = NA_REAL;
    return v;
")

# NaN, Inf, and NA all collapse into NA in the context of a logical vector
fun()
identical(fun(), c(FALSE, TRUE, FALSE, rep(NA, 3)))
