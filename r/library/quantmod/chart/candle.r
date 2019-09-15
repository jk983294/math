source("utils/utils.r")

dat <- read_stock_data()

chartSeries(dat, type = "candlesticks", theme = chartTheme("white"))
chartSeries(dat, type = "candlesticks", subset = "2018", theme = chartTheme("white"))
chartSeries(dat, type = "candlesticks", subset = "2018-05::2018-06", theme = chartTheme("white"))
