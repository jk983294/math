library(data.table)
library(glmnet)
library(tidyverse)

R2 <- function(pred, obs, formula = "corr", na.rm = FALSE) {
    n <- sum(complete.cases(pred))
    switch(formula, corr = cor(obs, pred, use = ifelse(na.rm, "complete.obs", "everything"))^2, traditional = 1 - (sum((obs - pred)^2, 
        na.rm = na.rm)/((n - 1) * var(obs, na.rm = na.rm))))
}

# dt1 <- read.csv('~/github/barn/train/factor.csv')
dt1 <- read.csv("~/github/MyTmp/sim/cmake-build-debug/pmfut/factor.csv")
plot(density(dt1$ret), main = "ret")
qqnorm(dt1$ret)
qqline(dt1$ret)

# split into training and validation set
split <- 0.8
rowTotal <- nrow(dt1)
rowSplit <- round(split * rowTotal)
index <- sample(rowTotal, rowSplit)
train <- dt1[index, ]
validation <- dt1[-index, ]

# fit model
fit <- lm(ret ~ . + 0, data = train)
summary(fit)

# out sample test
pred <- predict(fit, newdata = validation)
R2(pred, validation$ret)
pred_train <- predict(fit, newdata = train)
R2(pred_train, train$ret)
cor(cbind(pred, validation$ret))
cor(cbind(pred_train, train$ret))

# lasso
n.coef <- 100
lambdas <- 10^seq(-10, 10, length = n.coef)

fit <- glmnet(x, y, intercept = FALSE, lambda = lambdas)
coef.mat <- coef(fit)[-1, ]
oos.mat <- rep(0, n.coef)
for (i in 1:n.coef) {
    cur.coef <- coef.mat[, i]
    pred <- x_test %*% cur.coef
    oos.mat[i] <- R2(pred, y_test, "traditional")
}
plot(oos.mat, ylab = "R2", main = "ridge")

best <- which.max(oos.mat)
best
oos.mat[best]
coef.mat[, best]
