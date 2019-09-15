library(quantmod)
source("utils/utils.r")

dat <- read_stock_data()

WeekVolYa <- apply.weekly(Vo(dat), sum)  # sum from Monday to Friday
MonthVolYa <- apply.monthly(Vo(dat), sum)  # sum to month
QuarterVolYa <- apply.quarterly(Vo(dat), sum)  # sum to quarter
YearVolYa <- apply.yearly(Vo(dat), sum)  # sum to year
WeekAvgVolClYa <- apply.weekly(Vo(dat), mean)  # weekly average volume
