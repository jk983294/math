x <- rnorm(100)
y <- runif(100)

# method 1: check via histogram
hist(x, main = "hist", xlab = "Values", col = "blue", border = "red")
hist(y, main = "hist", xlab = "Values", col = "blue", border = "red")

# method 2: check qq-plot
qqnorm(x)
qqline(x, col = "red")

# method 3: Shapiro-Wilk test
# null hypothesis: the data is normally distributed
result <- shapiro.test(x)
result <- shapiro.test(y)

print(result)
