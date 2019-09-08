library(multcomp)
attach(cholesterol)
table(trt)
aggregate(response, by = list(trt), FUN = mean)
aggregate(response, by = list(trt), FUN = sd)
fit <- aov(response ~ trt)
summary(fit)
library(gplots)
plotmeans(response ~ trt, xlab = "Treatment", ylab = "Response", main = "Mean Plot\nwith 95% CI")
detach(cholesterol)

# Tukey HSD pairwise group comparisons
TukeyHSD(fit)
par(las = 2)
par(mar = c(5, 8, 4, 2))
plot(TukeyHSD(fit))
par(opar)

# Multiple comparisons
library(multcomp)
par(mar = c(5, 4, 6, 2))
tuk <- glht(fit, linfct = mcp(trt = "Tukey"))
plot(cld(tuk, level = 0.05), col = "lightgrey")
par(opar)
