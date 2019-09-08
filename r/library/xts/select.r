library(xts)

# prepare data
dat <- read.csv("../data/600848.csv", header = TRUE, row.names = "date")
dat <- xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))

dat["2019-03-07"]  # select by year month day
dat["2019-03"]  # select by year month
dat["2019"]  # select by year
dat["2019-03-07/2019-03-14"]  # select by range
dat["2019-03-07/2019-04"]  # select by range
dat[c("2019-03-07", "2019-03-11")]  # select by array

dat[endpoints(dat, on = "weeks", k = 2)]  # extract observations of every two week
dat[endpoints(dat, on = "years")]  # extract observations of every year

first(dat, "1 week")  # extract first 1 week
last(dat, "1 week")  # extract last 1 week
first(last(dat, "1 week"), "3 days")  # get first 3 days of the last week of data
