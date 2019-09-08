# Comparing nested models using the anova function
states <- as.data.frame(state.x77[, c("Murder", "Population", "Illiteracy", "Income", 
    "Frost")])
fit1 <- lm(Murder ~ Population + Illiteracy + Income + Frost, data = states)
fit2 <- lm(Murder ~ Population + Illiteracy, data = states)
anova(fit2, fit1)  # 检验不显著(p=0.994),不需要将这两个变量添加到线性模型中,可以将它们从模型中删除

# !!! ANOVA需要嵌套模型,而AIC方法不需要

# Comparing models with the AIC
fit1 <- lm(Murder ~ Population + Illiteracy + Income + Frost, data = states)
fit2 <- lm(Murder ~ Population + Illiteracy, data = states)
AIC(fit1, fit2)  # 优先选择AIC值较小的模型,它说明模型用较少的参数获得了足够的拟合度
