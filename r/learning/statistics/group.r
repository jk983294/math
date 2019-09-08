options(digits = 3)
attach(mtcars)
# group by (cyl, gear) and Calculate mean of each group for left columns
aggregate(mtcars, by = list(Group.cyl = cyl, Group.gear = gear), FUN = mean, na.rm = TRUE)

# by example
mystats <- function(x, na.omit = FALSE) {
    if (na.omit) 
        x <- x[!is.na(x)]
    m <- mean(x)
    n <- length(x)
    s <- sd(x)
    skew <- sum((x - m)^3/s^3)/n
    kurt <- sum((x - m)^4/s^4)/n - 3
    return(c(n = n, mean = m, stdev = s, skew = skew, kurtosis = kurt))
}
dstats <- function(x) sapply(x, mystats)
by(mtcars[c("mpg", "hp", "wt")], mtcars$am, dstats)  # group by am

# summaryBy
library(doBy)
summaryBy(mpg + hp + wt ~ am, data = mtcars, FUN = mystats)  # group by am

# describeBy
library(psych)
describeBy(mtcars[c("mpg", "hp", "wt")], list(am = mtcars$am)) # group by am
