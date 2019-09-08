# 通过一系列连续型或类别型预测变量来预测计数型结果变量

# look at dataset
data(breslow.dat, package = "robust")
names(breslow.dat)
summary(breslow.dat[c(6, 7, 8, 10)])

# plot distribution of post-treatment seizure counts
opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
attach(breslow.dat)
hist(sumY, breaks = 20, xlab = "Seizure Count", main = "Distribution of Seizures")
boxplot(sumY ~ Trt, xlab = "Treatment", main = "Group Comparisons")
par(opar)

# fit regression
fit <- glm(sumY ~ Base + Age + Trt, data = breslow.dat, family = poisson())
summary(fit)

# interpret model parameters
coef(fit)
exp(coef(fit))

# 当响应变量观测的方差比依据泊松分布预测的方差大时,泊松回归可能发生过度离势
# evaluate overdispersion
deviance(fit)/df.residual(fit)  # >> 1 means overdispersion
library(qcc)
qcc.overdispersion.test(breslow.dat$sumY, type = "poisson")  # 显著性检验的p<0.05,表明确实存在过度离势

# fit model with quasipoisson
fit.od <- glm(sumY ~ Base + Age + Trt, data = breslow.dat, family = quasipoisson())
summary(fit.od)
deviance(fit.od)/df.residual(fit.od)
