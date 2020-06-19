library(data.table)

dt <- data.table(school = c("NTU", "SMU", "NUS"), rank = c(2, 1, 3), size = c(1, 3, 2))
dim(dt)
class(dt$school)  # unlike data.frame, columns of character never converted to factors by default
