source("utils/utils.r")

stk <- read_stock_data()

periodicity(stk)  # Daily periodicity from 2017-03-07 to 2019-09-05
unclass(periodicity(stk))
to.weekly(stk)
to.monthly(stk)
periodicity(to.monthly(stk))
ndays(stk)  # 537
nweeks(stk)  # 113
nyears(stk)  # 3

getFX("USD/EUR")  # try some non-OHLC to start
periodicity(USDEUR)
to.weekly(USDEUR)
periodicity(to.weekly(USDEUR))
