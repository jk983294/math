# multi-variable statistics
x <- c(1, 3, 5, 2, 11, 9, 3, 9, 12, 3)
y <- c(4.4, 5.3, 7.2, 5.2, 8.5, 7.3, 6, 10.4, 10.2, 6.1)
cor(x, y)  # correlation 0.9075655
plot(x, y)

# covariance and correlation
states <- state.x77[, 1:6]
cov(states)  # 协方差
cor(states)  # pearson相关系数
cor(states, method = "spearman")  # spearman相关系数

# non square correlation
x <- states[, c("Population", "Income", "Illiteracy", "HS Grad")]
y <- states[, c("Life Exp", "Murder")]
cor(x, y)  # pearson相关系数

# 偏相关是指在控制一个或多个定量变量时,另外两个定量变量之间的相互关系
library(ggm)
# partial correlation of population and murder rate, controlling for income,
# illiteracy rate, and HS graduation rate
pcor(c(1, 5, 2, 3, 6), cov(states))
