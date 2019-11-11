now <- Sys.time()
typeof(now)  # 'double'
class(now)  # [1] 'POSIXct' 'POSIXt'
unclass(now)  # 1564303410 show its actual data

tomorrow <- unclass(now) + 24 * 60 * 60
class(tomorrow) <- c("POSIXct", "POSIXt")
tomorrow  # '2019-07-29 16:43:29 CST'

## POSIXct is the number of seconds since the epoch POSIXlt is a mixed text and character format
require(xts)
price_vector <- c(101.02, 101.03, 101.05)
dates <- c("03/12/2013 08:00:00.532123", "03/12/2013 08:00:05.700123", "03/12/2013 08:00:06.800123")
time_index <- strptime(dates, format = "%d/%m/%Y %H:%M:%OS")  # convert to POSIXlt
xts_price_vector <- xts(price_vector, time_index)

# time diff
time_vec <- as.POSIXct(dates, format = "%d/%m/%Y %H:%M:%OS")
diffs <- time_vec[-1] - time_vec[-length(time_vec)]
