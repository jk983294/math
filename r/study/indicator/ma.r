mySMA <- function(price, n) {
    sma <- c()
    sma[1:(n - 1)] <- NA
    for (i in n:length(price)) {
        sma[i] <- mean(price[(i - n + 1):i])
    }
    return(sma)
}
mySMA(1:10, 5)

myEMA <- function(price, n) {
    ema <- c()
    ema[1:(n - 1)] <- NA
    ema[n] <- mean(price[1:n])
    beta <- 2/(n + 1)
    for (i in (n + 1):length(price)) {
        ema[i] <- beta * price[i] + (1 - beta) * ema[i - 1]
    }
    return(ema)
}
myEMA(1:10, 5)

myRSI <- function(price, n) {
    U <- rep(0, length(price))
    D <- rep(0, length(price))
    rs <- rep(NA, length(price))
    rsi <- rep(NA, length(price))
    for (i in 2:length(price)) {
        if (price[i] > price[(i - 1)]) {
            U[i] <- 1
        } else if (price[i] < price[(i - 1)]) {
            D[i] <- 1
        }
        if (i > n) {
            if (sum(D[(i - n + 1):i]) == 0) {
                rsi[i] <- 100
            } else {
                rs[i] <- sum(U[(i - n + 1):i])/sum(D[(i - n + 1):i])
                rsi[i] <- rs[i]/(rs[i] + 1) * 100
            }
        }
    }
    return(rsi)
}
myRSI(sample(4:10, 20, replace = TRUE), 6)
