source("utils/utils.r")

dat <- read_stock_data()

# reChart takes most arguments of the original charting calls, and allows for quick modifications to your charts
candleChart(dat, subset = "2018", theme = chartTheme("white"))
reChart(major.ticks = "months", subset = "first 8 weeks")
