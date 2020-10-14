library(data.table)
library(glmnet)
library(tidyverse)

R2 <- function(pred, obs, formula = "corr", na.rm = FALSE) {
    n <- sum(complete.cases(pred))
    switch(formula, corr = cor(obs, pred, use = ifelse(na.rm, "complete.obs", "everything"))^2,
        traditional = 1 - (sum((obs - pred)^2, na.rm = na.rm)/((n - 1) * var(obs,
            na.rm = na.rm))))
}

split <- 0.8
dt1 <- read.csv("~/github/MyTmp/sim/cmake-build-debug/pmfut/factor.csv")
plot(density(dt1$ret), main = "ret")
qqnorm(dt1$ret)
qqline(dt1$ret)

# split into training and validation set
rowTotal = nrow(dt1)
rowSplit = round(split * rowTotal)
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

# ridge regression
y <- train$ret
x <- train %>% select(-ret) %>% data.matrix()
y_test <- validation$ret
x_test <- validation %>% select(-ret) %>% data.matrix()
n.coef <- 100
lambdas <- 10^seq(-4, 3, length = n.coef)

fit <- glmnet(x, y, intercept = FALSE, alpha = 0, lambda = lambdas)
summary(fit)
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

# ridge regression cross validation
y <- dt1 %>% select(ret) %>% scale(center = TRUE, scale = FALSE) %>% data.matrix()
X <- dt1 %>% select(-ret) %>% data.matrix()
fit <- cv.glmnet(X, y, intercept = FALSE, standardize=TRUE, alpha = 0, lambda = lambdas, nfolds = 10)
plot(fit)
lambda_best <- fit$lambda.min
model_cv <- glmnet(X, y, intercept = FALSE, standardize = TRUE, alpha = 0, lambda = lambda_best)
y_hat_cv <- predict(model_cv, X)
R2(y_hat_cv, y)
cor(cbind(y_hat_cv, y))

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
