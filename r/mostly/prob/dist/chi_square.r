# χ^2(n)
# χ^2(1) ~ Gamma(1/2, 1/2)
# χ^2(n) ~ Gamma(n/2, 1/2)
x <- (0L:50L) * 0.1
n <- 10
curve(dgamma(x, shape = n / 2, rate = 0.5), from = 0., to = 50)
