source("utils/utils.r")
library(TTR)

dat <- read_stock_data()

roc <- ROC(Cl(dat), n = 2)
head(roc, n = 5)

chartSeries(dat, type = "candlesticks", subset = "2018-05::2018-06", theme = chartTheme("white"))
addROC(n = 7)
