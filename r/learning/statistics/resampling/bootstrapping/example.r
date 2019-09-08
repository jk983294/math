# bootstrapping a single statistic
rsq <- function(formula, data, indices) {
    d <- data[indices, ]
    fit <- lm(formula, data = d)
    return(summary(fit)$r.square) # single statistic (R2)
}

library(boot)
set.seed(1234)
results <- boot(data = mtcars, statistic = rsq, R = 1000, formula = mpg ~ wt + disp)
print(results)
plot(results)
boot.ci(results, type = c("perc", "bca")) # 获取统计量的置信区间


# bootstrapping several statistics (regression coefficients)
bs <- function(formula, data, indices) {
    d <- data[indices, ]
    fit <- lm(formula, data = d)
    return(coef(fit))
}
library(boot)
set.seed(1234)
results <- boot(data = mtcars, statistic = bs, R = 1000, formula = mpg ~ wt + disp)
print(results)
# index1指截距项,index2指车重,index3指发动机排量
plot(results, index = 2)
boot.ci(results, type = "bca", index = 2)
boot.ci(results, type = "bca", index = 3)
