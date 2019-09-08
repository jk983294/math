library(xts)

# prepare data
dat <- read.csv("../data/600848.csv", header = TRUE, row.names = "date")
dat <- xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))

# remove duplicate
make.index.unique(dat, eps = 1e-04)
make.index.unique(dat, drop = TRUE)

# round index time
x <- Sys.time() + 1:1000
align.time(x, 10)  # every 10 seconds
align.time(x, 60)  # every minute
align.time(x, 10 * 60)  # every 10 min
