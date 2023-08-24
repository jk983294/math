library(data.table)
library(arrow)

dt <- read_feather("/tmp/20210915.feather")
dt1 <- dt[ins == 000021]
dt1[, .N]
write_feather(dt1, "/tmp/20210915.1.feather")