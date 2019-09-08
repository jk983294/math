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
plot(y, x, col = "blue", main = "Linear Regression", abline(lm(x ~ y)), cex = 1.3, 
    pch = 16)

# Nonlinear Least Square, training a model y = a1*x^2 + a2
x <- c(1.6, 2.1, 2, 2.23, 3.71, 3.25, 3.4, 3.86, 1.19, 2.21)
y <- c(5.19, 7.43, 6.94, 8.11, 18.75, 14.88, 16.06, 19.12, 3.21, 7.58)
plot(x, y)

init_param =  list(a1 = 1, a2 = 3)
model <- nls(y ~ a1 * x^2 + a2, start =init_param)  
print(model)
new_x <- data.frame(x = seq(min(x), max(x), len = 100))
lines(new_x$x, predict(model, new_x )) # plot fit curve

# get the sum of the squared residuals
print(sum(resid(model)^2))
# get the confidence intervals on the chosen values of the coefficients
print(confint(model))
