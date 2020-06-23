library(Rcpp)
library(inline)

Rcpp.package.skeleton("mypackage")
writeLines(system("tree", intern = TRUE))


rcpp_hello_world()
