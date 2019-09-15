source("utils/utils.r")
library(TTR)

dat <- read_stock_data()

M <- momentum(Cl(dat), n = 2)
head(M, n = 5)

chartSeries(dat, type = "candlesticks", subset = "2018-05::2018-06", theme = chartTheme("white"))
addMomentum(n = 1)
