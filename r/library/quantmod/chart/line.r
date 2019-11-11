source("utils/utils.r")

dat <- read_stock_data()

chartSeries(dat, type = "line", theme = chartTheme("white"))
chartSeries(dat, type = "line", subset = "2018", theme = chartTheme("white"))
chartSeries(dat, type = "line", subset = "2018-05::2018-06", theme = chartTheme("white"))

start_date <- "2018-05"
end_date <- "2018-06"
plot(dat[paste(start_date, "::", end_date, sep = "")])

# wrapper
lineChart(dat, subset = "2018", theme = chartTheme("white"))
