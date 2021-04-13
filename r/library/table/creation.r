library(data.table)

dt <- data.table(school = c("NTU", "SMU", "NUS"), rank = c(2, 1, 3), size = c(1, 3, 2))
dim(dt)
class(dt$school)  # unlike data.frame, columns of character never converted to factors by default

# add row
rbind(dt, list("a", 5, 6))

# empty dt
data <- data.table(va=numeric(), vb=numeric(), vc=numeric())

# create from list
l <- list(1:3, 4:6, 7:9)
dt <- setDT(lapply(l, unlist))