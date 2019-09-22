library(quantmod)
library(dplyr)

read_stock_data <- function(stock_path = "../data/600848.csv") {
    dat <- read.csv(stock_path, header = TRUE, row.names = "date")
    xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))
}

read_tushare_token <- function(path_ = "~/.ssh/tushare.token") {
    readLines(path_)
}

read_tushare_stock_daily <- function(stock_path = "~/data/index.399300.SZ.csv") {
    dat <- read.csv(stock_path, header = TRUE, row.names = "trade_date")
    dat <- dplyr::select(dat, -(X:ts_code))
    xts(dat, order.by = as.Date(rownames(dat), "%Y%m%d"))
}
