source("utils/utils.r")

stk <- read_stock_data()
idx <- read_index_data()

NS <- function(xdat) xdat/coredata(xdat)[1]
s <- NS(Cl(stk)) - 1
i <- NS(Cl(idx)) - 1
chartSeries(s, subset = "2018", theme = chartTheme("white"))
addTA(i, on = 1, col = "red", lty = "dotted")
