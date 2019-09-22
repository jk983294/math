source("utils/utils.r")

stk <- read_stock_data()

# general format for the above is CCYY-MM-DD HH:MM:SS
stk["2017"]  # returns all 2017 OHLC
stk["2018"]  # now just 2018
stk["2018-01"]  # now just January of 2018
stk["2017-06::2018-01-12"]  # range
stk["::"]  # everything in stk
stk["2018::"]  # everything in stk, from 2018 onward
non.contiguous <- c("2019-01", "2019-02", "2019-12")
stk[non.contiguous]

last(stk)  # returns the last observation.
last(stk, 8)  # returns the last 8 observation
# let's try something a bit cooler.
last(stk, "3 weeks")
last(stk, "-3 weeks")  # all except the last 3 weeks
last(stk, "3 months")
last(first(stk, "2 weeks"), "3 days")
