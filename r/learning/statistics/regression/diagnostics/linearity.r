states <- as.data.frame(state.x77[, c("Murder", "Population", "Illiteracy", "Income", 
    "Frost")])
fit <- lm(Murder ~ Population + Illiteracy + Income + Frost, data = states)

# Assessing linearity
library(car)
crPlots(fit)  # 若图形存在非线性,则说明你可能对预测变量的函数形式建模不够充分, 那么就需要添加一些曲线成分

# Global test of linear model assumptions
library(gvlma)
gvmodel <- gvlma(fit)
summary(gvmodel)
