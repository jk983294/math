source("utils/utils.r")
require(urca)

price_distribution <- function(symbol) {
    dat <- read_stock_data(symbol)
    prices <- dat$close
    mean_prices <- round(mean(prices), 2)
    sd_prices <- round(sd(prices), 2)
    hist(prices, breaks = 100, prob = T, cex.main = 0.9)
    abline(v = mean_prices, lwd = 2)
    legend("topright", cex = 0.8, border = NULL, bty = "n", paste("mean=", mean_prices, "; sd=", sd_prices))
}

return_distribution <- function(symbol) {
    dat <- read_stock_data(symbol)
    prices <- dat$close
    returns <- diff(log(prices))  # compute log returns
    hist(returns, breaks = 100, prob = T, cex.main = 0.9)  # returns is more stationary than prices
}

# result shows return is more normally distributed than price
price_distribution("../data/000001.csv")
return_distribution("../data/000001.csv")

check_price_stationary <- function(symbol) {
    dat <- read_stock_data(symbol)
    prices <- dat$close
    test <- ur.kpss(as.numeric(prices))  # null hypothesis that the time series is stationary
    test@teststat  # test statistic > critical value, we reject null hypothesis
    test@cval  # critical values
    test
}

check_return_stationary <- function(symbol) {
    dat <- read_stock_data(symbol)
    prices <- dat$close
    returns <- diff(log(prices))
    test_returns <- ur.kpss(as.numeric(returns))
    test_returns@teststat
    test_returns@cval
    test_returns
}

# result shows return is more stationary than price
check_price_stationary("../data/000001.csv")
check_return_stationary("../data/000001.csv")
