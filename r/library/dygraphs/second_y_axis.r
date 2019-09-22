source("utils/utils.r")
library(dplyr)

dat <- read.csv("../data/600848.csv", header = TRUE, row.names = "date")
dat <- dplyr::select(dat, open, volume)
stk <- xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))

dygraph(stk, main = "data") %>% dySeries("volume", axis = "y2") %>% dyOptions(colors = RColorBrewer::brewer.pal(2, "Set2")) %>% dyHighlight(highlightSeriesOpts = list(strokeWidth = 3)) %>% 
    dyRangeSelector(height = 20)
