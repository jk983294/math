library(Rcpp)

Rcpp.package.skeleton(name = "cpppkg", cpp_files = c("cpppkg.cpp"))

pkgbuild::compile_dll()  # precompile, run befor devtools::document()
devtools::document() # Update the NAMESPACE and doc

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

MyTimesTwo(42)
ben = new(Student, name = "Ben", age = 26, male = TRUE)
ben$LikesBlue()
ben$GetNumbers()

sam = new(Teacher, name = "sam", age = 26, male = TRUE)
sam$LikesBlue()

jk = new(BadStudent, name = "jk", age = 26, male = TRUE)
jk$GetNumbers()
#jk$LikesBlue()
