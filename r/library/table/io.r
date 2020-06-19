library(data.table)

dt <- data.table(x = c("a", "b", "c", "b"), y = c(2, 1, 3, 4), z = c(1, 3, 2, 2))
fwrite(dt, "/tmp/test_dt.csv")

dt1 <- fread("/tmp/test_dt.csv")
