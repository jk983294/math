source("utils/utils.r")
library(TTR)

dat <- read_stock_data()

chartSeries(dat, theme = "white", TA = NULL)

# one way
chartSeries(dat, theme = "white", TA = "addVo();addBBands();addCCI()")

# another way
chartSeries(dat, theme = "white")  #draw the chart
addVo()  # add volume
addBBands()  # add Bollinger Bands
addCCI()  # add Commodity Channel Index
