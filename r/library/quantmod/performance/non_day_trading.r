library(PerformanceAnalytics)
source("utils/utils.r")

stk <- read_stock_data()

qty <- 300
day <- 14

signal <- c()  # trade signal
signal[1:(day + 1)] <- 0

price <- Cl(stk)

stock <- c()  # stock holding
stock[1:(day + 1)] <- 0

cash <- c()
cash[1:(day + 1)] <- 10000

# generate trading signal
rsi <- RSI(price, day)  #rsi is the lag of RSI
for (i in (day + 1):length(price)) {
    if (rsi[i] < 30) {
        # buy one more unit if rsi < 30
        signal[i] <- 1
    } else if (rsi[i] < 50) {
        # no change if rsi < 50
        signal[i] <- 0
    } else {
        # sell if rsi > 50
        signal[i] <- -1
    }
}
signal <- reclass(signal, price)

# assume buying at closing price. We keep track of how cash and stock changes
trade <- Lag(signal)  # rsi is the lag of RSI
for (i in (day + 1):length(price)) {
    if (trade[i] >= 0) {
        stock[i] <- stock[i - 1] + qty * trade[i]
        cash[i] <- cash[i - 1] - qty * trade[i] * price[i]
    } else {
        stock[i] <- 0
        cash[i] <- cash[i - 1] + stock[i - 1] * price[i]
    }
}
stock <- reclass(stock, price)
cash <- reclass(cash, price)

# To evaluate performance, we calculate equity using cash and stock holdings.
equity <- c()
equity[1:(day + 1)] <- 10000

return <- c()
return[1:(day + 1)] <- 0

for (i in (day + 1):length(price)) {
    equity[i] <- stock[i] * price[i] + cash[i]
    return[i] <- equity[i]/equity[i - 1] - 1
}
equity <- reclass(equity, price)
return <- reclass(return, price)
charts.PerformanceSummary(return, main = "Non-Day-Trading")

chart_Series(equity, main="equity line")
chart_Series(cash, name="Cash Holding") # check the cash account over time
chart_Series(stock, name="Stock Holding") 
