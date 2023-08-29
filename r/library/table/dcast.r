library(data.table)

II <- c(1,2,1,2)
Time <- c(1,1,2,2)
signal1 <- c(1.1,2.2,3.3,4.4)
signal2 <- c(1.2,2.3,3.4,4.5)
signals <- data.table(Time,II,signal1, signal2)
setkey(signals, Time, II)

# convert to long form
(ldt <- melt(signals, id.vars = c("Time", "II"), measure.vars = c("signal1", "signal2")))

# long form to wide form
(wdt <- dcast(ldt, Time + II ~ variable, value.var = c("value")))

# check equal with original table
all.equal(signals, wdt, ignore.row.order = TRUE, ignore.col.order = TRUE)
