# A random vector follows a multivariate normal distribution if every linear combination
# of its components has a univariate normal distribution.
# X=(X_1, ..., X_n)^T
# pdf: f(x)= 1 / ((2π)^n/2 * ∣Σ∣^1/2) exp(-1/2 * (x−μ)^T * Σ^−1 * (x−μ))

library(mvtnorm)

mean_vec <- c(0, 0)
rho <- 0.7 # correlation
cov_mat <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
r <- rmvnorm(n = 10^3, mean = mean_vec, sigma = cov_mat) # draw sample
plot(r)
cov(r)

gen_mvnorm <- function(n, rho) {
    tau <- sqrt(1 - rho^2)
    x <- rnorm(n)
    y <- rnorm(n)
    z <- x
    w <- rho * x + tau * y
    cbind(z, w)
}
r <- gen_mvnorm(10^3, rho)
cov(r)
