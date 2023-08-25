x <- c("why", "video", "cross", "extra", "deal", "authority")
y <- c("a", "b", "c", "d", "e", "f")

str_c(x, collapse = ", ") # concat one list
str_c(x, y, sep = "-") # zip
str_c(x, y, sep = "-", collapse = ", ") # zip then collapse to one line
str_flatten(x, collapse = ",") # collapse into a single string
str_dup(x, times = 2L) # every word get repeated