# Listing 8.15 - Function for k-fold cross-validated R-square
shrinkage <- function(fit, k = 10) {
    require(bootstrap)
    
    # define functions
    theta.fit <- function(x, y) {
        lsfit(x, y)
    }
    theta.predict <- function(fit, x) {
        cbind(1, x) %*% fit$coef
    }
    
    # matrix of predictors
    x <- fit$model[, 2:ncol(fit$model)]
    # vector of predicted values
    y <- fit$model[, 1]
    
    results <- crossval(x, y, theta.fit, theta.predict, ngroup = k)
    r2 <- cor(y, fit$fitted.values)^2  # raw R2
    r2cv <- cor(y, results$cv.fit)^2  # cross-validated R2
    cat("Original R-square =", r2, "\n")
    cat(k, "Fold Cross-Validated R-square =", r2cv, "\n")
    cat("Change =", r2 - r2cv, "\n")
}

# using it
states <- as.data.frame(state.x77[, c("Murder", "Population", "Illiteracy", "Income", 
    "Frost")])
fit <- lm(Murder ~ Population + Income + Illiteracy + Frost, data = states)
shrinkage(fit)
fit2 <- lm(Murder ~ Population + Illiteracy, data = states)
shrinkage(fit2)
