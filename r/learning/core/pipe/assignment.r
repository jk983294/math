library(magrittr)
x <- c(0.109, 0.359, 0.63, 0.996, 0.515, 0.142, 0.017, 0.829, 0.907)

# compound assignment operator %<>%
x %<>% abs %>% sort  # x <- sort(abs(x))

# %<>% needs to be the first pipe operator in the chain for this to work
