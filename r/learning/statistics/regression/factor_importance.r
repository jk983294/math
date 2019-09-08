# Calculating standardized regression coefficients
# 根据预测变量与响应变量的相关系数来进行排序
states <- as.data.frame(state.x77[, c("Murder", "Population", "Illiteracy", "Income", 
    "Frost")])
zstates <- as.data.frame(scale(states))  # 数据标准化为均值为0、标准差为1的数据集
zfit <- lm(Murder ~ Population + Income + Illiteracy + Frost, data = zstates)
coef(zfit)

# relative weight for calculating relative importance of predictors
relweights <- function(fit, ...) {
    R <- cor(fit$model)
    nvar <- ncol(R)
    rxx <- R[2:nvar, 2:nvar]
    rxy <- R[2:nvar, 1]
    svd <- eigen(rxx)
    evec <- svd$vectors
    ev <- svd$values
    delta <- diag(sqrt(ev))
    lambda <- evec %*% delta %*% t(evec)
    lambdasq <- lambda^2
    beta <- solve(lambda) %*% rxy
    rsquare <- colSums(beta^2)
    rawwgt <- lambdasq %*% beta^2
    import <- (rawwgt/rsquare) * 100
    import <- as.data.frame(import)
    row.names(import) <- names(fit$model[2:nvar])
    names(import) <- "Weights"
    import <- import[order(import), 1, drop = FALSE]
    dotchart(import$Weights, labels = row.names(import), xlab = "% of R-Square", 
        pch = 19, main = "Relative Importance of Predictor Variables", sub = paste("Total R-Square=", 
            round(rsquare, digits = 3)), ...)
    return(import)
}

# Applying the relative weight function
states <- as.data.frame(state.x77[, c("Murder", "Population", "Illiteracy", "Income", 
    "Frost")])
fit <- lm(Murder ~ Population + Illiteracy + Income + Frost, data = states)
relweights(fit, col = "blue")
