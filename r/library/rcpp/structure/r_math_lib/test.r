# PKG_CXXFLAGS=`Rscript -e 'Rcpp:::CxxFlags()'` PKG_LIBS=`Rscript -e
# 'Rcpp:::LdFlags()'` R CMD SHLIB mynorm.cpp rm -f mynorm.o mynorm.so

dyn.load("./library/rcpp/structure/r_math_lib/mynorm.so")
.Call("mypnorm", 1:4)
