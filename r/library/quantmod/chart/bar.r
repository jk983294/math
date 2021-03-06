source("utils/utils.r")

dat <- read_stock_data()

chartSeries(dat, type = "bar", theme = chartTheme("white"))
chartSeries(dat, type = "bar", subset = "2018", theme = chartTheme("white"))
chartSeries(dat, type = "bar", subset = "2018-05::2018-06", theme = chartTheme("white"))

# wrapper
barChart(dat, subset = "2018", theme = chartTheme("white"))
