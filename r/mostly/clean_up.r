# remove all variable
rm(list = ls())

# handle NA vec
(x <- c(1L, 2L, NA, 4L, NA, 5L))
(bad <- is.na(x))
(x[!bad])

# handle NA multi-vec
(y <- c("a", "b", NA, "d", NA, "f"))
(good <- complete.cases(list(x, y)))
(x[good])
(y[good])

# clear packages
detach("package:datasets", unload = TRUE)

# clear console
cat("\014") # ctrl+L