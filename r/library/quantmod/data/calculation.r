source("utils/utils.r")

stk <- read_stock_data()

# changes for all combinations of columns
open_close_change <- OpCl(stk)  # (close(i) - open(i)) / open(i)
open_open_change <- OpOp(stk)  # (open(i) - open(i - 1)) / open(i - 1)
high_close_change <- HiCl(stk)  # (close(i) - high(i)) / high(i)

lag1 <- Lag(Cl(stk))  # x[i] = x[i-1]
lag135 <- Lag(Cl(stk), c(1, 3, 5))  # One, three, and five period lags
next1 <- Next(Cl(stk))  # x[i] = x[x+1]
delta1 <- Delt(Op(stk), Cl(stk), k = 1)  # (close[x+k] - open[x]) / open[x]
delta123 <- Delt(Op(stk), Cl(stk), k = 1:3)  # Open to close one-day, two-day and three-day lags

# calculate periodic minimums, maximums, sums, and products
apply.weekly(stk, FUN = function(x) {
    max(Cl(x))
})
period.apply(stk, endpoints(stk, on = "weeks"), FUN = function(x) {
    max(Cl(x))
})
period.max(Cl(stk), endpoints(stk, on = "weeks"))  # same thing - only 50x faster
# period.min, period.sum, period.prod, period.max

# Period Returns
dr <- dailyReturn(stk)  # returns by day, (close[i] - open[i]) / open[i]
wr <- weeklyReturn(stk)  # returns by week, (close[last day of week] - open[first day of week]) / open[first day of week]
mr <- monthlyReturn(stk)  # returns by month, indexed by yearmon
qr <- quarterlyReturn(stk)
yr <- yearlyReturn(stk)
allReturns(stk)  # daily,weekly,monthly,quarterly, and yearly
