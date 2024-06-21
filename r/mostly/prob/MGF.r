# MGF of N(u, sigma^2)
M <- function(t, mean = 0, sd = 1) {
    exp(mean * t + sd^2 * t^2 / 2.)
}
curve(M, from = -3, to = 3)

# continuous moment using numerical integration
# 2th moments of N(0, 1)
g <- function(x) x^2 * dnorm(x)
integrate(g, lower = -Inf, upper = Inf)

# discrete moment using sum
# 2th moments of pois(7)
g <- function(k) k^2 * dpois(k, 7)
sum(g(0:100))
