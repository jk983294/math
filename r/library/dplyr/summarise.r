# summarise() reduces multiple values down to a single summary
library(dplyr)

df <- data.frame(ID = c("a", "b", "a", "b", "a", "b"), x = c(1, 2, 3, 4, 5, 6), y = c(1, 4, 2, 3, 6, 2))
summarize(df, count = n(), xbar = mean(x), yvar = var(y))

# measure of location
summarize(df, count = n(), ymean = mean(y), ymedian = median(y))

# measure of spread
summarize(df, count = n(), ysd = sd(y), yIQR = IQR(y), ymad = mad(y))

# measure of rank
summarize(df, count = n(), ymin = min(y), yquantile25 = quantile(y, 0.25), ymax = max(y))

# measure of position
summarize(df, count = n(), yfirst = first(y), y2 = nth(y, 2), ylast = last(y))

# count
summarize(df, count = n(), non_na_count = sum(!is.na(x)), unique_count = n_distinct(x))
summarize(df, condition_count = sum(x < 3))
