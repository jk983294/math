source("utils/utils.r")
library(dplyr)

dat <- read.csv("../data/600848.csv", header = TRUE, row.names = "date")
dat <- dplyr::select(dat, open, close, high, low)
stk <- xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))

dygraph(stk, main = "data") %>%
    dyOptions(colors = RColorBrewer::brewer.pal(4, "Set2")) %>% 
    dyHighlight(highlightSeriesOpts = list(strokeWidth = 3)) %>% 
    dyRangeSelector(height = 20)

# dySeries give specific option for one line
dygraph(stk, main = "data") %>%
    dyOptions(colors = RColorBrewer::brewer.pal(4, "Set2")) %>% 
    dySeries("low", stepPlot = TRUE, fillGraph = TRUE, color = "red") %>%
    dyRangeSelector(height = 20)
