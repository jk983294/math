# Permutation tests for simple linear regression
library(lmPerm)
set.seed(1234)
fit <- lmp(weight ~ height, data = women, perm = "Prob")
summary(fit)


# Permutation tests for polynomial regression
library(lmPerm)
set.seed(1234)
fit <- lmp(weight ~ height + I(height^2), data = women, perm = "Prob")
summary(fit)

# Permutation tests for multiple regression
library(lmPerm)
set.seed(1234)
states <- as.data.frame(state.x77)
fit <- lmp(Murder ~ Population + Illiteracy + Income + Frost, data = states, perm = "Prob")
summary(fit)
