library(quantmod)

getSymbols("AAPL", src = "yahoo")
getSymbols("AAPL", from = "2012-01-01", to = "2019-08-25")
head(AAPL)

AAPL <- last(AAPL, "1 year")
AAPL <- first(AAPL, "2 months")

Open <- Op(AAPL)  # Open Price
High <- Hi(AAPL)  # High price
Low <- Lo(AAPL)  # Low price
Close <- Cl(AAPL)  # Close Price
Volume <- Vo(AAPL)  # Volume
AdjClose <- Ad(AAPL)  # Adjusted close

dat <- read.csv("../data/600848.csv", header = TRUE, row.names = "date")
dat <- xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))
dat <- last(dat, "1 year")
Open <- Op(dat)  # Open Price
High <- Hi(dat)  # High price
Low <- Lo(dat)  # Low price
Close <- Cl(dat)  # Close Price
Volume <- Vo(dat)  # Volume
