library(dplyr)

(df <- data.frame(ID = c("a", "b", "a", "b", "a", "b", "a", "b"), x = c(1, 2, 3, 4, 5, 6, 7, 8), y = c(1, 4, 2, 3, 6, 2, 8, -1)))
(newdf <- group_by(df, ID))
summarize(newdf, count = n(), xbar = mean(x), yvar = var(y))

# create id number by group
(df <- df %>% group_by(ID) %>% dplyr::mutate(id_num = cur_group_id()))

# calculate mean by group and add to df
(df <- df %>% dplyr::group_by(ID) %>% dplyr::mutate(gr_mean_x = mean(x)) %>% as.data.frame())

# calculate cumsum by group and add to df
(df <- df %>% dplyr::group_by(ID) %>% dplyr::mutate(cs = cumsum(x)) %>% as.data.frame())

# center
(ans <- df %>% group_by(ID) %>% dplyr::summarise(mean = mean(x)))
(ans <- df %>% group_by(ID) %>% dplyr::summarise(median = median(x)))

# spread
(ans <- df %>% group_by(ID) %>% dplyr::summarise(sd = sd(x)))
(ans <- df %>% group_by(ID) %>% dplyr::summarise(IQR = IQR(x)))
(ans <- df %>% group_by(ID) %>% dplyr::summarise(mad = mad(x)))

# range
(ans <- df %>% group_by(ID) %>% dplyr::summarise(min = min(x)))
(ans <- df %>% group_by(ID) %>% dplyr::summarise(max = max(x)))
(ans <- df %>% group_by(ID) %>% dplyr::summarise(quantile = quantile(x, 0.25)))

# top n
(ans <- df %>% arrange(desc(y)) %>% group_by(ID) %>% slice(1:2))

# max/min
(ans <- df %>% group_by(ID) %>% slice_max(n = 1, y))
(ans <- df %>% group_by(ID) %>% slice_min(n = 1, y))

# position
(ans <- df %>% group_by(ID) %>% dplyr::summarise(first = first(x)))
(ans <- df %>% group_by(ID) %>% dplyr::summarise(last = last(x)))
(ans <- df %>% group_by(ID) %>% dplyr::summarise(nth = nth(x, 2)))

# count
(ans <- df %>% group_by(ID) %>% dplyr::summarise(n = n()))
(ans <- df %>% group_by(ID) %>% dplyr::summarise(n_distinct = n_distinct(x)))

# logical
(ans <- df %>% group_by(ID) %>% dplyr::summarise(val = any(x > 2)))
(ans <- df %>% group_by(ID) %>% dplyr::summarise(val = all(x > 2)))