source("utils/utils.r")

dat <- read_stock_data()

candleChart(dat, , theme = chartTheme("white", up.col = "red", dn.col = "green"))
zoomChart("last 2 week")  # zoom in
reChart(subset = "last 2 months")  # zoom out
zoomChart("last 1 months")  # zoom in
zoomChart()  # zoom out
zooom()  # interactive example
