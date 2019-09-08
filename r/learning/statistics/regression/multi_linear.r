# y = a + b1x1 + b2x2 +...bnxn
df <- data.frame(matrix(rnorm(18, 1), ncol = 3))
model <- lm(X3 ~ X1 + X2, data = df)
print(model)

# use model to predict
a <- data.frame(matrix(rnorm(2, 1), ncol = 2))
result <- predict(model, a)
print(result)

# examining bivariate relationships
states <- as.data.frame(state.x77[, c("Murder", "Population", "Illiteracy", "Income", 
    "Frost")])
cor(states)
library(car)
scatterplotMatrix(states, spread = FALSE, smoother.args = list(lty = 2), main = "Scatter Plot Matrix")


# Multiple linear regression
fit <- lm(Murder ~ Population + Illiteracy + Income + Frost, data = states)
summary(fit)
confint(fit) # 95%的置信区间

# Mutiple linear regression with interaction term
fit <- lm(mpg ~ hp + wt + hp:wt, data = mtcars)
summary(fit)

library(effects)
plot(effect("hp:wt", fit, , list(wt = c(2.2, 3.2, 4.2))), multiline = TRUE)
