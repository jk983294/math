library(Rcpp)

code <- "
    #include <gsl/gsl_const_mksa.h>
    // decl of constants

    std::vector<double> volumes() {
        std::vector<double> v(5);
        v[0] = GSL_CONST_MKSA_US_GALLON; // 1 US gallon
        v[1] = GSL_CONST_MKSA_CANADIAN_GALLON; // 1 Canadian gallon
        v[2] = GSL_CONST_MKSA_UK_GALLON; // 1 UK gallon
        v[3] = GSL_CONST_MKSA_QUART; // 1 quart
        v[4] = GSL_CONST_MKSA_PINT; // 1 pint
        return v;
    }"

gslVolumes <- cppFunction(code, depends = "RcppGSL")
gslVolumes()
