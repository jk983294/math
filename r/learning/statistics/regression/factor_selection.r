# 逐步回归中,模型会一次添加或者删除一个变量,直到达到某个判停准则为止

# Backward stepwise selection 一次删除一个变量直到会降低模型质量为止
library(MASS)
states <- as.data.frame(state.x77[, c("Murder", "Population", "Illiteracy", "Income", 
    "Frost")])
fit <- lm(Murder ~ Population + Illiteracy + Income + Frost, data = states)
stepAIC(fit, direction = "backward")  # smaller AIC is better


# All subsets regression, 全子集回归是指所有可能的模型都会被检验
library(leaps)
states <- as.data.frame(state.x77[, c("Murder", "Population", "Illiteracy", "Income", 
    "Frost")])
leaps <- regsubsets(Murder ~ Population + Illiteracy + Income + Frost, data = states, 
    nbest = 4)
plot(leaps, scale = "adjr2")  # 通过R平方、调整R平方或Mallows Cp统计量等准则来选择“最佳”模型
library(car)
subsets(leaps, statistic = "cp", main = "Cp Plot for All Subsets Regression")
abline(1, 1, lty = 2, col = "red")  # 越好的模型离截距项和斜率均为1的直线越近
