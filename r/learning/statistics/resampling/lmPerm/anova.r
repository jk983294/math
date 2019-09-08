# Permutation test for One-Way ANOVA
library(lmPerm)
library(multcomp)
set.seed(1234)
fit <- lmp(response ~ trt, data = cholesterol, perm = "Prob")
anova(fit)


# Permutation test for One-Way ANCOVA
library(lmPerm)
set.seed(1234)
fit <- lmp(weight ~ gesttime + dose, data = litter, perm = "Prob")
anova(fit)


# Permutation test for Two-way ANOVA
library(lmPerm)
set.seed(1234)
fit <- lmp(len ~ supp * dose, data = ToothGrowth, perm = "Prob")
anova(fit)
