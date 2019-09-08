library(xts)

# prepare data
dat <- read.csv("../data/600848.csv", header = TRUE, row.names = "date")
dat <- xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))

dat["2017-03-07"] <- 4  # replace by index value
dat[2] <- 2  # replace by index
