library(purrr)

y <- list(1, 2, 3)
z <- list(4, 5, 6)
l1 <- list(x = "a", y = "z")
l2 <- list(x = "a", y = "z")

purrr::map2_dbl(y, z, ~ .x / .y)
purrr::map2_dbl(y, z, `+`)
purrr::map2_chr(l1, l2, paste, collapse = ",", sep = ":")
purrr::map2_lgl(l2, l1, `%in%`)
purrr::walk2(objs, paths, save) # return invisibly

purrr::map2(y, z, function(a, b) a / b)
