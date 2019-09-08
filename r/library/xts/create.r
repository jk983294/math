library(xts)

data <- rnorm(5)
dates1 <- seq(as.Date("2017-05-01"), length = 5, by = "days")
dates2 <- Sys.Date() - 1:10
dates3 <- as.POSIXct(Sys.Date() + 1:10)
xts1 <- xts(x = 1:10, order.by = dates2)
xts2 <- xts(x = data, order.by = dates1)
xts3 <- xts(x = rnorm(10), order.by = dates3, born = as.POSIXct("1899-05-08"))
xts4 <- xts(x = 1:10, order.by = Sys.Date() + 1:10)

# convert from df
data(AirPassengers)
xts5 <- as.xts(AirPassengers)

# copy the timestamp from xts2 to xts6
xts6 <- c(2, 1, 0, -1, -2)
xts6 <- reclass(xts6, xts2)
xts6
