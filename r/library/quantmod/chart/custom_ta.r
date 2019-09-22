source("utils/utils.r")
library(TTR)

dat <- read_stock_data()

chartSeries(dat, type = "candlesticks", subset = "2018-05::2018-06", theme = chartTheme("white"))
sma <- SMA(Cl(dat), n = 14)
addTA(sma, on = 1, col = "red")

# Using newTA it is possible to create your own generic TA function
addOpCl <- newTA(OpCl, col = "green", type = "h")
addOpCl()
