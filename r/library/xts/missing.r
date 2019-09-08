library(xts)

# prepare data
dat <- read.csv("../data/600848.csv", header = TRUE, row.names = "date")
dat <- xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))

dat[2] <- NA

na.omit(dat)
dat_prev <- na.locf(dat)  # Fill missing values using prev observation
dat_next <- na.locf(dat, fromLast = TRUE)  # Fill missing values using next observation
dat_approx <- na.approx(dat)  # Interpolate NAs
