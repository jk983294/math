source("utils/utils.r")
library(TTR)

dat <- read_stock_data()

rsi = RSI(Cl(dat), n = 14)
tail(rsi, n = 5)

chartSeries(dat, type = "candlesticks", subset = "2018-05::2018-06", theme = chartTheme("white"))
addRSI(n = 14, maType = "EMA")
