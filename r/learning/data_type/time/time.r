now <- Sys.time()
typeof(now)  # 'double'
class(now)  # [1] 'POSIXct' 'POSIXt'
unclass(now)  # 1564303410 show its actual data

tomorrow <- unclass(now) + 24 * 60 * 60
class(tomorrow) <- c("POSIXct", "POSIXt")
tomorrow  # '2019-07-29 16:43:29 CST'
