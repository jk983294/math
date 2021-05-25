# filter() picks cases based on their values
library(dplyr)

df <- data.frame(Price = c(1.2, 2.4, 3.6, 4.8), month = c(1, 2, 1, 2), day = c(1, 2, 3, 4), desc = c("a", "b", "c", "d"))
dplyr::filter(df, month == 1, day == 3)
dplyr::filter(df, month == 1 | day == 2)
dplyr::filter(df, desc %in% c("a", "b"))
dplyr::filter(df, !(Price > 3))

# sample data
df <- data_frame(a = c(1, 2, 3, NA), b = c(5, Inf, 8, 8), c = c(9, 10, Inf, 11), d = c('a', 'b', 'c', 'd'))

# across all columns
df %>% filter_all(all_vars(!is.infinite(.)))

# note that is.finite() does not work with NA or strings
df %>% filter_all(all_vars(is.finite(.)))

# checking only numeric columns
df %>% filter_if(~is.numeric(.), all_vars(!is.infinite(.)))

# checking only select columns, in this case a through c
df %>% filter_at(vars(a:c), all_vars(!is.infinite(.)))
