# rm(list = ls())  # 删除其他函数，保证package只包含下面两个函数

# myfun <- function(x) { x + 7 }

# myfun2 <- function(x) { x * 7 }

# 运行下面命令会将上面两个函数（已运行）导出到package
package.skeleton(name = "mypkg")

# check with load_all

library(devtools)
library(roxygen2)
my.pkg <- as.package("mypkg")
load_all(my.pkg)
devtools::document(my.pkg)

# build in terminal
$ R CMD build mypkg
$ R CMD INSTALL mypkg_1.0.tar.gz
