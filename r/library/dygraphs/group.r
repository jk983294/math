source("utils/utils.r")
library(dplyr)

dat <- read.csv("../data/600848.csv", header = TRUE, row.names = "date")
open_ <- dplyr::select(dat, open)
volume_ <- dplyr::select(dat, volume)
stk_open <- xts(open_, order.by = as.Date(rownames(open_), "%Y-%m-%d"))
stk_volume <- xts(volume_, order.by = as.Date(rownames(volume_), "%Y-%m-%d"))

dygraph(stk_open, main = "open", group = "group1")
dygraph(stk_volume, main = "volume", group = "group1")
