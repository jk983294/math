library(data.table)

II <- c(1,2,1,2)
Time <- c(1,1,2,2)
signal1 <- c(1.1,2.2,3.3,4.4)
signal2 <- c(1.2,2.3,3.4,4.5)
signals <- data.table(Time,II,signal1, signal2)
dcast(signals, Time ~ II, value.var=list("signal1"))
