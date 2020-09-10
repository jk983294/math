library(Rcpp)
library(RcppArmadillo)

RcppArmadillo.package.skeleton("cpppkg")

# create project on existing folder cpppkg, run below in R
compileAttributes(verbose=TRUE)
# src/RcppExports.cpp updated
# R/RcppExports.R updated

package_native_routine_registration_skeleton(dir="/home/kun/github/math/r/cpppkg/")
# copy result into src/init.c

# compile package
$ R CMD build cpppkg
$ R CMD INSTALL cpppkg_1.0.tar.gz

# test
library(cpppkg)

ben = new(Student, name = "Ben", age = 26, male = TRUE)
ben$LikesBlue()