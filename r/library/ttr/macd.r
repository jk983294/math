source("utils/utils.r")
library(TTR)

dat <- read_stock_data()

macd <- MACD(Cl(dat), nFast = 12, nSlow = 26, nSig = 9, maType = EMA)
tail(macd, n = 5)

chartSeries(dat, type = "candlesticks", subset = "2018-05::2018-06", theme = chartTheme("white"))
addMACD(fast = 12, slow = 26, signal = 9, type = "EMA")
