dt <- data.table(x = c("a", "b", "a", "b"), y = c(2., 1., 3., NA), z = c(1., 3., NA, 2.), m = c("m", "m", "n", "m"))

my_subset <- function(data, col, val) {
    data[data[[col]] == val, ]
}
my_subset1 <- function(data, col, val) {
    # deparse(substitute(...)) to catch the actual names of objects passed to function
    data <- deparse(substitute(data))
    col  <- deparse(substitute(col))
    val  <- paste0("'", val, "'")
    text <- paste0("subset(", data, ", ", col, " == ", val, ")")
    # subset(data, col == val)
    eval(parse(text = text)[[1L]])
}

my_subset2 <- function(data, col, val) {
  eval(substitute(subset(data, col == val)))
}

my_subset(dt, col = "x", val = "a")
my_subset1(dt, x, "a")
my_subset2(dt, x, "a")


square <- function(x) x * x
dt[, outer(inner(var1) + inner(var2)),
    env = list(outer = "sqrt", inner = "square", var1 = "y", var2 = "z")]
dt[filter_col %in% filter_val,
   .(var1, var2, out = outer(inner(var1) + inner(var2))),
   by = by_col,
   env = list(
     outer = "sqrt",
     inner = "square",
     var1 = "y",
     var2 = "z",
     out = "calc_val",
     filter_col = "x",
     filter_val = I(c("a", "b")),
     by_col =  "m"
  )]
