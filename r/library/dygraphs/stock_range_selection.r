source("utils/utils.r")

stk <- read_stock_data()

price <- OHLC(dat)
dygraph(price) %>% dyCandlestick() %>% dyRangeSelector(height = 20)
