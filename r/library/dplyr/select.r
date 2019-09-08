# select() picks variables based on their names
library(dplyr)

df <- data.frame(x = c("NTU", "SMU", "NUS"), rank = c(2, 1, 3), size = c(1, 3, 2))
dplyr::select(df, x, rank)

# select all columns between x and rank (inclusive)
select(df, x:rank)

# select all columns except those from year to day (inclusive)
select(df, -(x:rank))

select(df, size, everything())  # move size to the first column
