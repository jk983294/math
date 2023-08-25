library(purrr)

# https://github.com/rstudio/cheatsheets/blob/main/purrr.pdf

x <- list(a = 1L:10L, b = 11L:20L, c = 21L:30L)
l1 <- list(x = c("a", "b"), y = c("c", "d"))
purrr::map(l1, sort, decreasing = TRUE)
purrr::map_dbl(x, mean)  # return a double vector
purrr::map_int(x, length) # return a int vector
purrr::map_chr(l1, paste, collapse = "") # return a char vector
purrr::map_lgl(x, is.integer) # return a logic vector
purrr::walk(x, print)  # return invisibly

purrr::map(x, function(x) sum(x))
purrr::modify(x, function(x) x + 1.)
purrr::modify_if(x, is.numeric, function(x) x + 1.)
purrr::reduce(x, sum)
purrr::compact(x)  # rm empty elements
purrr::keep(x, is.numeric)
purrr::every(x, is.numeric)
purrr::some(x, is.numeric)
purrr::none(x, is.numeric)
purrr::has_element(x$a, 5L)
purrr::pluck(x$a, 5L)  # extract
purrr::assign_in(x$a, 5L, -5L)
purrr::modify_in(x$a, 5L, function(x) x - 5L)