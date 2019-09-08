# 离群点是指那些模型预测效果不佳的观测点

states <- as.data.frame(state.x77[, c("Murder", "Population", "Illiteracy", "Income", 
    "Frost")])
fit <- lm(Murder ~ Population + Illiteracy + Income + Frost, data = states)

# Assessing outliers
library(car)
outlierTest(fit)

# 高杠杆值观测点,即与其他预测变量有关的离群点。换句话说,它们是由许多异常的预测变量值组合起来的,与响应变量值没有关系.
# 帽子均值为p/n,其中p是模型估计的参数数目(包含截距项),n是样本量。一般来说,若观测点的帽子值大于帽子均值的2或3倍,就可以认定为高杠杆值点.
# Identifying high leverage points
hat.plot <- function(fit) {
    p <- length(coefficients(fit))
    n <- length(fitted(fit))
    plot(hatvalues(fit), main = "Index Plot of Hat Values")
    abline(h = c(2, 3) * p/n, col = "red", lty = 2)
    identify(1:n, hatvalues(fit), names(hatvalues(fit)))
}
hat.plot(fit)  # press ESC to quit

# 强影响点,即对模型参数估计值影响有些比例失衡的点.
cutoff <- 4/(nrow(states) - length(fit$coefficients) - 2)  # Cook距离
plot(fit, which = 4, cook.levels = cutoff)
abline(h = cutoff, lty = 2, col = "red")

# 变量添加图,即对于每个预测变量X_k,绘制X_k在其他k–1个预测变量上回归的残差值相对于响应变量在其他k–1个预测变量上回归的残差值的关系图.
# Added variable plots add id.method='identify' to interactively identify points
library(car)
avPlots(fit, ask = FALSE, id.method = "identify")

# Influence Plot 将离群点、杠杆值和强影响点的信息整合到一幅图形
library(car)
influencePlot(fit, id.method = "identify", main = "Influence Plot", sub = "Circle size is proportial to Cook's Distance")
