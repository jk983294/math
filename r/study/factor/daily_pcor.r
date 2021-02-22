library(data.table)
library(ggplot2)

dt1 <- fread("/tmp/dump_pcor.csv")
dim(dt1)

ggplot(dt1, aes(di)) + geom_line(aes(y = f7, colour = "y"))
