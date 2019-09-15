source("utils/utils.r")
library(TTR)

dat <- read_stock_data()

sma <- SMA(Cl(dat), n = 20)
head(sma, n = 5)
tail(sma, n = 5)

chartSeries(dat, type = "candlesticks", subset = "2018-05::2018-06", theme = chartTheme("white"))
addSMA(n = 30, on = 1, col = "blue")
addSMA(n = 200, on = 1, col = "red")

ema <- EMA(Cl(dat), n = 20)
head(sma, n = 5)
tail(ema, n = 5)

chartSeries(dat, type = "candlesticks", subset = "2018-05::2018-06", theme = chartTheme("white"))
addEMA(n = 30, on = 1, col = "blue")
addEMA(n = 200, on = 1, col = "red")
