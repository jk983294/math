# X ~ binp(n, p), if n*p= const, use 泊松分布
n <- 100000L
p <- 1e-5
lambda <- n * p
pbinom(1L, size = n, prob = p) # check two cdfs
ppois(1L, lambda = lambda)

# X ~ binp(n, p), if n*p*(1-p) > 10, use 正态分布
# p越接近0.5,近似效果越好,因为对称性, 形态越像正态
# p越接近两端, 形态越skew, 需要越大的n使得接近正态
n <- 1000000L
p <- 2e-4
(sd_ <- sqrt(n * p * (1. - p)))
(mean_ <- n * p)
x <- 210L
pbinom(x, size = n, prob = p) # check two cdfs
pnorm(x, mean = mean_, sd = sd_)
