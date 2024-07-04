nsim <- 1e4
n <- 12L
x <- matrix(runif(n * nsim), nrow = nsim, ncol = n) # converge fast
x <- matrix(rexp(n * nsim, rate = 1.), nrow = nsim, ncol = n) # converge slow
xbar <- rowMeans(x)
hist(xbar)
