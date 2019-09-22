library(PerformanceAnalytics)
source("utils/utils.r")

stk <- read_stock_data()

# Buy signal based on EMA rule. Sell signal based on RSI rule. Tie-breaking: buy-signal has priority
n <- 14
delta <- 0.005
price <- Cl(stk)
r <- price/Lag(stk) - 1
rsi <- RSI(price, n)
signal <- c()  # first signal is NA
signal[1:n] <- 0


# Generate Trading Signal
for (i in (n + 1):length(price)) {
    if (r[i] > delta) {
        signal[i] <- 1
    } else if (rsi[i] > 70) {
        signal[i] <- -1
    } else signal[i] <- 0
}
signal <- reclass(signal, price)

## Apply Trading Rule
trade3 <- Lag(signal)
ret3 <- dailyReturn(stk) * trade3
names(ret3) <- "Combine"
charts.PerformanceSummary(ret3, main = "RSI")
