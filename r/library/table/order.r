library(dplyr)
library(data.table)

data <- data.table(x1 = c("a", "b", "c", "b"), x2 = c(2, 1, 3, 4), x3 = c(1, 3, 2, 2))

data[order(x2), ]
(setorder(data, x1, -x2))
(setkey(data, x1, x2))
data[order(-x2), ] # reverse order
dplyr::arrange(data, x2) # order data with dplyr

# multi dimension sort
data[order(x1, -x3)]

# set secondary index for further calculation
(setindex(data, x2))