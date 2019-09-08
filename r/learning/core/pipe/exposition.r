library(magrittr)
# exposition pipe operator %$%
data.frame(z = rnorm(100)) %$% ts.plot(z)  # %$% make sure that data frame are exposed to ts.plot, similar to with statement
