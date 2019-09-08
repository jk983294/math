library(xts)

# prepare data
dat <- read.csv("../data/600848.csv", header = TRUE, row.names = "date")
dat <- xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))

# data attributes
time(dat[1])  # get time stamp
as.numeric(dat[1])  # get the value
core_data <- coredata(dat)  # 'matrix'
index_data <- index(dat)  # 'Date'
indexClass(dat)  # 'Date'
indexClass(convertIndex(dat, "POSIXct"))  # replacing index class
indexFormat(dat) <- "%Y/%m/%d"
.index(dat)  # raw numeric index, second from epoch
time(dat)  # get time stamp
.indexwday(dat)  # value of weekday
.indexhour(dat)  # value of hour
start(dat)  # first time
end(dat)  # last time

# timezone
indexTZ(dat)  # 'UTC'
tzone(dat) <- "Asia/Hong_Kong"
tzone(dat)  # 'Asia/Hong_Kong'

# periods
periodicity(dat)  # Daily periodicity from 2017-03-07 to 2019-09-05
to.yearly(dat)  # yearly OHLC
to.quarterly(dat)  # quarterly OHLC
to.monthly(dat)  # monthly OHLC
to.weekly(dat)  # weekly OHLC
to.period(dat, period = "quarters")  # quarterly OHLC
to.period(dat, period = "years")  # yearly OHLC
nweeks(dat)  # count the weeks
nmonths(dat)  # count the months
nquarters(dat)  # count the quarters
nyears(dat)  # count the years
