library(quantmod)

read_stock_data <- function(stock_path = "../data/600848.csv") {
    dat <- read.csv(stock_path, header = TRUE, row.names = "date")
    xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))
}
