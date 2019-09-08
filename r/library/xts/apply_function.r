library(xts)

# prepare data
dat <- read.csv("../data/600848.csv", header = TRUE, row.names = "date")
dat <- xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))

# choose observations
ep1 <- endpoints(dat, on = "weeks", k = 2)  # index of every two week
ep2 <- endpoints(dat, on = "years")  # index of year
dat[ep1]  # extract observations of every two week
dat[ep2]  # extract observations of every year
dat_yearly <- split(dat, f = "years")  # split into several lists

# apply function
period.apply(dat, INDEX = ep2, FUN = mean)  # yearly mean
lapply(dat_yearly, FUN = mean)  # yearly mean

# find the last 1 month observation in each year
do.call(rbind, lapply(split(dat, "years"), function(w) last(w, n = "1 month")))

# calculate cumulative sum
do.call(rbind, lapply(split(dat, "years"), cumsum))

# apply standard deviation to rolling window
rollapply(dat, 3, sd)

# Arithmetic Operations
xts1 <- xts(x = rnorm(5), order.by = Sys.Date() + (1:5 + 2))
xts2 <- xts(x = 1:5, order.by = Sys.Date() + 1:5)
xts1 + as.numeric(xts2)  # Addition
xts1 * as.numeric(xts2)  # Multiplication
coredata(xts1) - xts2  # Subtraction
coredata(xts1)/xts2  # Division

# lag
lag(xts1, 2)
xts1 - lag(xts1)  # period-over-period differences, xts1[i] - xts1[i- 1]
diff(xts1, lag = 1, differences = 1)  # lagged differences

# Reindexing
xts1 + merge(xts2, index(xts1), fill = 100)  # Addition
xts1 - merge(xts2, index(xts1), fill = na.locf)  # Subtraction
## na.locf Last Observation Carried Forward

# Merging
merge(xts2, xts1, join = "inner")  # inner join of xts2 and xts1
merge(xts2, xts1, join = "left", fill = 0)  # left join of xts2 and xts1,
rbind(xts1, xts2)
cbind(xts1, xts2)
