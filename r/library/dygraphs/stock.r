library(quantmod)
library(dygraphs)

dat <- read.csv("../data/600848.csv", header = TRUE, row.names = "date")
dat <- xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))
price <- OHLC(dat)
head(price, n = 3)

graph <- dygraph(price)

# Shading
dyShading(graph, from = "2018-04-09", to = "2018-10-11", color = "#FFE6E6")

# Event line
graph <- dyEvent(graph, "2017-6-29", "a event", labelLoc = "bottom")
dyEvent(graph, "2019-5-6", "b event", labelLoc = "bottom")

# Candle Chart
selected_dat <- tail(dat, n = 30)
graph <- dygraph(OHLC(selected_dat))
dyCandlestick(graph)
