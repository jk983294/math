x <- seq(0, 50, by = 1)
pmf <- dbinom(x, 50, 0.5)
plot(x, pmf)

# P(X <= 25) with n = 50, p = 0.5
prob <- pbinom(25, 50, 0.5) # 0.5561376
# P(X <= a) = prob
quantile <- qbinom(prob, 50, 0.5) # 25

# sample based on binomial distribution
samples <- rbinom(1000, 50, 0.5)
hist(samples)
