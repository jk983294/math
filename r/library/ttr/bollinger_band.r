source("utils/utils.r")
library(TTR)

dat <- read_stock_data()

bb <- BBands(Cl(dat), s.d = 2)
tail(bb, n = 5)

chartSeries(dat, type = "candlesticks", subset = "2018-05::2018-06", theme = chartTheme("white"))
addBBands(n = 20, sd = 2)
