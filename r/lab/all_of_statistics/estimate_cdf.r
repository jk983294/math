x <- seq(-3, 3, by = 0.05)
cdf <- pnorm(x, mean = 0, sd = 1)  # Cumulative Distribution Function
plot(x, cdf)

# create a sample of 50 numbers which are normally distributed
sample <- rnorm(100)
hist(sample)
Fn <- ecdf(sample)  # Empirical Cumulative Distribution Function
plot(Fn)
ecdf_vals <- Fn(x)
plot(ecdf_vals)

# confidence band
n <- length(x)
alpha <- 0.05
epsilon_n <- sqrt(1/(2 * n) * log(2/alpha))
conf_lower_bound <- pmax(ecdf_vals - epsilon_n, 0)
conf_upper_bound <- pmin(ecdf_vals + epsilon_n, 1)

plot(ecdf_vals)
lines(conf_lower_bound)
lines(conf_upper_bound)
lines(cdf)
