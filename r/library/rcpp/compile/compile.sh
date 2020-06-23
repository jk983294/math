PKG_CXXFLAGS="-I/usr/local/lib/R/site-library/Rcpp/include" PKG_LIBS="-L/usr/local/lib/R/site-library/Rcpp/libs -lRcpp" R CMD SHLIB fibonacci.cpp


PKG_CXXFLAGS=`Rscript -e 'Rcpp:::CxxFlags()'` PKG_LIBS=`Rscript -e 'Rcpp:::LdFlags()'` R CMD SHLIB fibonacci.cpp

rm -f fibonacci.o fibonacci.so
