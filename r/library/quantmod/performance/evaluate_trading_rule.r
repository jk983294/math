library(PerformanceAnalytics)
source("utils/utils.r")

stk <- read_stock_data()

# simple filter rule suggests buying when the price increases a lot compared to the yesterday price
price <- Cl(stk)  # close price
r <- price/Lag(price) - 1  # % price change
delta <- 0.005  #threshold
signal <- c(0)  # first date has no signal

for (i in 2:length(price)) {
    # Loop over all trading days (except the first)
    if (r[i] > delta) {
        signal[i] <- 1
    } else if (r[i] < -delta) {
        signal[i] <- -1
    } else signal[i] <- 0
}

signal <- reclass(signal, price)  # assign time to action variable using reclass

# Charting with Trading rule
chartSeries(stk, type = "line", theme = chartTheme("white"))
addTA(signal, type = "S", col = "red")

trade <- Lag(signal, 1)  # trade based on yesterday signal
ret <- dailyReturn(stk) * trade  # buy at open, sell at close, trading size: all in

names(ret) <- "Naive"
charts.PerformanceSummary(ret, main = "Naive")

# rsi based
day <- 14
signal <- c()  # initialize vector
rsi <- RSI(price, day)  # rsi is the lag of RSI
signal[1:day + 1] <- 0  # 0 because no signal until day+1

for (i in (day + 1):length(price)) {
    if (rsi[i] < 30) {
        # buy if rsi < 30
        signal[i] <- 1
    } else {
        # no trade all if rsi > 30
        signal[i] <- 0
    }
}
signal <- reclass(signal, Cl(price))
trade2 <- Lag(signal)

# construct a new variable ret2
ret2 <- dailyReturn(stk) * trade2
names(ret2) <- "RSI"
charts.PerformanceSummary(ret2, main = "RSI")

# Now compare strategies :
retall <- cbind(ret, ret2)
charts.PerformanceSummary(retall, main = "Naive v.s. RSI")
