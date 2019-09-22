source("utils/utils.r")

stk <- read_stock_data()

Op(stk)  # extract Open
Hi(stk)  # extract High
Lo(stk)  # extract Low
Cl(stk)  # extract
Vo(stk)  # extract Volume

is.OHLC(stk)  # check if is OHLC
has.OHLC(stk)  # check if has OHLC
has.Op(stk)  # check if has Open
has.Cl(stk)  # check if has Close
has.Hi(stk)  # check if has High
has.Lo(stk)  # check if has Low
has.Vo(stk)  # check if has Volume

seriesHi(stk)  # where and what was the high point 
seriesLo(stk) # where and what was the low point
seriesHi(Cl(stk))
seriesLo(Cl(stk))
