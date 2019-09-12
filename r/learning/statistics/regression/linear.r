# Linear Regression y = ax + b
x <- c(151, 174, 138, 186, 128, 136, 179, 163, 152, 131)
y <- c(63, 81, 56, 91, 47, 57, 76, 72, 62, 48)
model <- lm(y ~ x)
print(model)
print(summary(model))  # error, residual etc

# use model to predict
a <- data.frame(x = 170)
result <- predict(model, a)
print(result)

# visualize
plot(y, x, col = "blue", main = "Linear Regression", abline(lm(x ~ y)), cex = 1.3, pch = 16)

# another example
fit <- lm(weight ~ height, data = women)
summary(fit)  # R-squared:  0.991 模型可以解释体重99.1%的方差
women$weight
fitted(fit)  # use model recalculated weight, compare against women$weight
residuals(fit)  # error, women$weight - fitted(fit)
plot(women$height, women$weight, main = "Women Age 30-39", xlab = "Height (in inches)", ylab = "Weight (in pounds)")
# add the line of best fit
abline(fit)

library(car)
scatterplot(weight ~ height, data = women, spread = FALSE, smoother.args = list(lty = 2), pch = 19, main = "Women Age 30-39", xlab = "Height (inches)", 
    ylab = "Weight (lbs.)")

# lmtest
library(sandwich)
library(lmtest)

x <- c(1, 2, 3, 5, 6, 7, 10, 12, 13)
y <- c(1, 4, 5, 6, 7, 8, 9, 10, 15)
z <- c(2, 3, 7, 8, 9, 12, 8, 7, 6)
df <- data.frame(x = x, y = y, z = z)
MLR <- lm(y ~ x + z, df)
coeftest(MLR, vcov = vcovHC(MLR, "HC1"))
