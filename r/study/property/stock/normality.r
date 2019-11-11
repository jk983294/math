source("utils/utils.r")

## superimpose such a normal distribution onto our empirical daily return data
dat <- read_stock_data()
prices <- dat$close
returns <- diff(log(prices))  # compute log returns
mu <- mean(returns, na.rm = TRUE)
sigma <- sd(returns, na.rm = TRUE)
x <- seq(-5 * sigma, 5 * sigma, length = nrow(returns))
hist(returns, breaks = 100, main = "Histogram of returns", cex.main = 0.8, prob = TRUE)
lines(x, dnorm(x, mu, sigma), col = "red", lwd = 2)  # leptokyrtic distribution

# QQ plot for normality check
par(mfrow = c(1, 2))
qqnorm(as.numeric(returns), main = "empirical returns qqplot()", cex.main = 0.8)
qqline(as.numeric(returns), lwd = 2)
grid()
normal_data <- rnorm(nrow(returns), mean = mu, sd = sigma)  # Normal random data
qqnorm(normal_data, main = "Normal returns", cex.main = 0.8)
qqline(normal_data, lwd = 2)
grid()

# other statistical tests: Kolmogorov-Smirnov, Shapiro-Wilks, Cramer-von Mises, Anderson-Darling
answer <- shapiro.test(as.numeric(returns))
answer
