library(data.table)
library(glmnet)
library(tidyverse)
library(doParallel)
registerDoParallel(4)

dt1 <- read.csv("~/github/barn/train/factor.csv")

# lasso regression
y <- train$ret
x <- train %>% select(-ret) %>% data.matrix()
y_test <- validation$ret
x_test <- validation %>% select(-ret) %>% data.matrix()
n.coef <- 100
lambdas <- 10^seq(-4, 3, length = n.coef)

fit <- glmnet(x, y, intercept = FALSE, alpha = 1, lambda = lambdas)
summary(fit)
coef.mat <- coef(fit)[-1, ]
oos.mat <- rep(0, n.coef)
for (i in 1:n.coef) {
    cur.coef <- coef.mat[, i]
    pred <- x_test %*% cur.coef
    oos.mat[i] <- R2(pred, y_test, "traditional")
}
plot(oos.mat, ylab = "R2", main = "lasso")

best <- which.max(oos.mat)
best
oos.mat[best]

# lasso regression cross validation
y <- dt1 %>% select(ret) %>% scale(center = TRUE, scale = FALSE) %>% data.matrix()
X <- dt1 %>% select(-ret) %>% data.matrix()
fit <- cv.glmnet(X, y, intercept = FALSE, standardize = TRUE, alpha = 1, lambda = lambdas, nfolds = 10, parallel = TRUE)
plot(fit)
lambda_best <- fit$lambda.min
model_cv <- glmnet(X, y, intercept = FALSE, standardize = TRUE, alpha = 1, lambda = lambda_best)
y_hat_cv <- predict(model_cv, X)
R2(y_hat_cv, y)
cor(cbind(y_hat_cv, y))

# See how increasing lambda shrinks the coefficients -------------------------- Each line shows coefficients for one variables, for
# different lambdas.  The higher the lambda, the more the coefficients are shrinked towards zero.
res <- glmnet(X, y, alpha = 1, lambda = lambdas, standardize = FALSE)
plot(res, xvar = "lambda")
legend("bottomright", lwd = 1, col = 1:6, legend = colnames(X), cex = 0.7)
