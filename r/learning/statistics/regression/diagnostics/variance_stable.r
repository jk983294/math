states <- as.data.frame(state.x77[, c("Murder", "Population", "Illiteracy", "Income", 
    "Frost")])
fit <- lm(Murder ~ Population + Illiteracy + Income + Frost, data = states)

# Assessing homoscedasticity 判断误差方差是否恒定
library(car)
ncvTest(fit)  # 检验不显著(p=0.19),说明满足方差不变假设
spreadLevelPlot(fit)  # Suggested power transformation close to 1 means 方差恒定
