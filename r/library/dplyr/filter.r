# filter() picks cases based on their values
library(dplyr)

df <- data.frame(Price = c(1.2, 2.4, 3.6, 4.8), month = c(1, 2, 1, 2), day = c(1, 2, 3, 4), desc = c("a", "b", "c", "d"))
dplyr::filter(df, month == 1, day == 3)
dplyr::filter(df, month == 1 | day == 2)
dplyr::filter(df, desc %in% c("a", "b"))
dplyr::filter(df, !(Price > 3))
