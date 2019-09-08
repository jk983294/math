# arrange() changes the ordering of the rows
library(dplyr)

df <- data.frame(school = c("NTU", "SMU", "NUS"), rank = c(2, 1, 3), size = c(1, 3, 2))
dplyr::arrange(df, rank, size)
dplyr::arrange(df, desc(rank), size)
