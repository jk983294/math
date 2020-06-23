fun <- cxxfunction(signature(), plugin = "Rcpp", body = "
    Rcpp::CharacterVector v(3);
    v[0] = \"The quick brown\";
    v[1] = \"fox\";
    v[2] = R_NaString;
    return v;
")

# Character vectors can also be converted to std::vector<std::string>
fun()
