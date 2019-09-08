library(foreach)
library(doParallel)
registerDoParallel(cores = 4)

eig <- function(n, p) {
    x <- matrix(rnorm(1e+05), ncol = 100)
    r <- cor(x)
    eigen(r)$values
}
n <- 1e+06
p <- 100
k <- 500

# %do% 操作符按顺序运行函数, 函数使用 %dopar% 操作符进行并行运算
system.time(res <- foreach(i = 1:k, .combine = rbind) %do% eig(n, p))
system.time(res <- foreach(i = 1:k, .combine = rbind) %dopar% eig(n, p))
# .combine=rbind操作符追加对象res作为行
