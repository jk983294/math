# study the effect of the categorical variable by using it along with the
# predictor variable and comparing the regression lines for each level of the
# categorical variable. Such an analysis is termed as Analysis of Covariance also
# called as ANOVA
y = floor(runif(9, min = 0, max = 2))  # categorical variable {0, 1}
df <- data.frame(cbind(matrix(rnorm(18, 1), ncol = 2), y))
colnames(df) <- c("X1", "X2", "y")
model1 <- aov(y ~ X1 + X2, data = df)
print(summary(model1))
model2 <- aov(y ~ X1 * X2, data = df)
print(summary(model2))

print(anova(model1, model2))  # compare the two models
