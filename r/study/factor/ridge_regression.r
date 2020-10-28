library(data.table)
library(glmnet)
library(tidyverse)
library(doParallel)
registerDoParallel(4)

dt1 <- read.csv("~/github/barn/train/factor.csv")

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
fit <- cv.glmnet(X, y, intercept = FALSE, standardize = TRUE, alpha = 0, lambda = lambdas, nfolds = 10, parallel = TRUE)
plot(fit)
lambda_best <- fit$lambda.min
model_cv <- glmnet(X, y, intercept = FALSE, standardize = TRUE, alpha = 0, lambda = lambda_best)
y_hat_cv <- predict(model_cv, X)
R2(y_hat_cv, y)
cor(cbind(y_hat_cv, y))

# Use information criteria to select lambda
X_scaled <- scale(X)
n.coef <- 10
lambdas <- 10^seq(-4, 3, length = n.coef)
foreach(index = seq(lambdas)) %dopar% {
    # Run model
    model <- glmnet(X, y, alpha = 0, lambda = lambdas[index], standardize = TRUE)
    # Extract coefficients and residuals (remove first row for the intercept)
    betas <- as.vector((as.matrix(coef(model))[-1, ]))
    resid <- y - (X_scaled %*% betas)
    # Compute hat-matrix and degrees of freedom
    ld <- lambdas[index] * diag(ncol(X_scaled))
    H <- X_scaled %*% solve(t(X_scaled) %*% X_scaled + ld) %*% t(X_scaled)
    df <- tr(H)
    # Compute information criteria
    my_ret_list <- list()
    my_ret_list$aic <- nrow(X_scaled) * log(t(resid) %*% resid) + 2 * df
    my_ret_list$bic <- nrow(X_scaled) * log(t(resid) %*% resid) + 2 * df * log(nrow(X_scaled))
    return(my_ret_list)
}

aic <- sapply(result, function(x) x$aic)
bic <- sapply(result, function(x) x$bic)

# Plot information criteria against tried values of lambdas
plot(log(lambdas), aic, col = "orange", type = "l", ylim = c(190, 260), ylab = "Information Criterion")
lines(log(lambdas), bic, col = "skyblue3")
legend("bottomright", lwd = 1, col = c("orange", "skyblue3"), legend = c("AIC", "BIC"))
