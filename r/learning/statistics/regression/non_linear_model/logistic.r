# 通过一系列连续型和/或类别型预测变量来预测二值型结果变量

# get summary statistics
data(Affairs, package = "AER")
summary(Affairs)
table(Affairs$affairs)

# convert to binary variable
Affairs$ynaffair[Affairs$affairs > 0] <- 1
Affairs$ynaffair[Affairs$affairs == 0] <- 0
Affairs$ynaffair <- factor(Affairs$ynaffair, levels = c(0, 1), labels = c("No", "Yes"))
table(Affairs$ynaffair)

# fit full model
fit.full <- glm(ynaffair ~ gender + age + yearsmarried + children + religiousness + 
    education + occupation + rating, data = Affairs, family = binomial())
summary(fit.full)  # big p value, less significant, since null hypothesis is coefficient = 0

# fit reduced model based on p vaule of fit.full
fit.reduced <- glm(ynaffair ~ age + yearsmarried + religiousness + rating, data = Affairs, 
    family = binomial())
summary(fit.reduced)

# diagnostics
plot(predict(fit.full, type = "response"), residuals(fit.full, type = "deviance"))
plot(hatvalues(fit.full))
plot(rstudent(fit.full))
plot(cooks.distance(fit.full))

plot(predict(fit.reduced, type = "response"), residuals(fit.reduced, type = "deviance"))
plot(hatvalues(fit.reduced))
plot(rstudent(fit.reduced))
plot(cooks.distance(fit.reduced))

# compare models
anova(fit.reduced, fit.full, test = "Chisq")  # 卡方值不显著(p=0.21),表明四个预测变量的新模型与九个完整预测变量的模型拟合程度一样好

# interpret coefficients 解释模型参数
coef(fit.reduced)
exp(coef(fit.reduced))

# calculate probability of extramariatal affair by marital ratings
testdata <- data.frame(rating = c(1, 2, 3, 4, 5), age = mean(Affairs$age), yearsmarried = mean(Affairs$yearsmarried), 
    religiousness = mean(Affairs$religiousness))
testdata$prob <- predict(fit.reduced, newdata = testdata, type = "response")
testdata

# calculate probabilites of extramariatal affair by age
testdata <- data.frame(rating = mean(Affairs$rating), age = seq(17, 57, 10), yearsmarried = mean(Affairs$yearsmarried), 
    religiousness = mean(Affairs$religiousness))
testdata$prob <- predict(fit.reduced, newdata = testdata, type = "response")
testdata

# 当出现过度离势时,仍可使用glm函数拟合Logistic回归,但此时需要将二项分布改为类二项分布

# evaluate overdispersion
overdispersion <- deviance(fit.reduced)/df.residual(fit.reduced)

# compare two model to check overdispersion
fit <- glm(ynaffair ~ age + yearsmarried + religiousness + rating, family = binomial(), 
    data = Affairs)
fit.od <- glm(ynaffair ~ age + yearsmarried + religiousness + rating, family = quasibinomial(), 
    data = Affairs)
pchisq(summary(fit.od)$dispersion * fit$df.residual, fit$df.residual, lower = F)
