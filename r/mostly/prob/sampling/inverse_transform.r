## Universality of the Uniform

# simulate exponential r.v. using uniform r.v.
# Uniform U, can be plugged to inverse of cumulative density function and we would have r.v. X, in accordance to the CDF.
# U ~ Unif(0, 1), F is CDF
# let X = F^-1(U), then X ~ F
# P(X ≤ x) = P(F^-1(U) ≤ x) = P(F(F^-1(U)) ≤ F(x)) = P(U ≤ F(x)) = F(x)
n <- 10000L
lambda <- 2.
u <- runif(n, min = 0, max = 1)
plot(density(u))
FM::stats(u)
inv_exp_cdf <- function(v, lambda) {
    -log(1. - v) / lambda
}
exp_rv <- inv_exp_cdf(u, 2.)
FM::stats(exp_rv) # check mean = 1/lambda, sd = 1/lambda
plot(density(exp_rv))
hist(exp_rv)
