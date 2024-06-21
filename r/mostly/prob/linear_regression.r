num <- 100L
x <- rnorm(num)
y <- 3 + 5 * x + rnorm(num)

b <- cov(x, y) / var(x)
a <- mean(y) - b * mean(x)

plot(x, y)
abline(a = a, b = b) # y = a + b * x
