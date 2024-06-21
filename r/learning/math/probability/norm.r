# Normal Distribution
x <- seq(-5, 5, by = 0.1)
density <- dnorm(x, mean = 0, sd = 1)
plot(x, density)
curve(dnorm, from=-3, to=3, n=1000)

cdf <- pnorm(x, mean = 0, sd = 1)  # Cumulative Distribution Function
plot(x, cdf)

# quantile whose cumulative value matches the probability value
x <- seq(0, 1, by = 0.02)  # probability range from [0, 1]
quantile <- qnorm(x, mean = 0, sd = 1)
plot(x, quantile)

# create a sample of 50 numbers which are normally distributed
sample <- rnorm(1000)
hist(sample, main = "Normal DIstribution")
