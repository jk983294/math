library(PerformanceAnalytics)
source("utils/utils.r")

stk <- read_stock_data()

qty <- 1000  # Trade unit: 1000 stocks per trade
day <- 14

signal <- c()  # trade signal with size
signal[1:(day + 1)] <- 0

price <- Cl(stk)

wealth <- c()
wealth[1:(day + 1)] <- 1e+06  # Wealth: 1 million

return <- c()
return[1:(day + 1)] <- 0

profit <- c()
profit[1:(day + 1)] <- 0

# generate trading signal with size
rsi <- RSI(price, day)  # rsi is the lag of RSI
for (i in (day + 1):length(price)) {
    if (rsi[i] < 30) {
        # buy one more unit if rsi < 30
        signal[i] <- signal[i - 1] + 1
    } else if (rsi[i] < 50) {
        # no change if rsi < 50
        signal[i] <- signal[i - 1]
    } else {
        # sell if rsi > 50
        signal[i] <- 0
    }
}
signal <- reclass(signal, price)

# apply Trade Rule
Close <- Cl(stk)
Open <- Op(stk)
trade <- Lag(signal)
for (i in (day + 1):length(price)) {
    profit[i] <- qty * trade[i] * (Close[i] - Open[i])
    wealth[i] <- wealth[i - 1] + profit[i]
    return[i] <- (wealth[i]/wealth[i - 1]) - 1
}
ret3 <- reclass(return, price)

charts.PerformanceSummary(ret3, main = "Trade Size")
