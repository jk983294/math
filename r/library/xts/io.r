library(xts)

file_name <- "../data/600848.csv"

# read file
dat <- read.csv(file_name, header = TRUE, row.names = "date")
class(dat)  # 'data.frame'
dat <- xts(dat, order.by = as.Date(rownames(dat), "%Y-%m-%d"))
class(dat)  # 'xts' 'zoo'

# write zoo
write.zoo(dat, index.name = "date", sep = ",", file = "/tmp/600848.csv")

# read zoo
dat2 <- read.csv.zoo("/tmp/600848.csv", index.column = list("date"), sep = ",")
