# common paradigm may be to compute and create the core C++ aspects of an object
# at the C++ level and to complement the objects at the R language level

f1 <- cxxfunction(signature(x = "any"), plugin = "Rcpp", body = "
    RObject y(x) ;
    List res(3) ;
    res[0] = y.isS4();
    res[1] = y.hasSlot(\"z\");
    res[2] = y.slot(\"z\");
    return res;
")

f2 <- cxxfunction(signature(x = "any"), plugin = "Rcpp", body = "
    S4 foo(x);
    foo.slot(\".Data\") = \"foooo\";
    return foo;
")
