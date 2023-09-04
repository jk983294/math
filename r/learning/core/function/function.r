# function and invoke
f <- function(x = 5) {
    x + 1
}
f(2)  # by position
f(x = 2)  # by name
f()  # default value

args(f)  # argument string
formals(f)  # argument list
body(f)

roll <- function(bones = 1:6) {
    dice <- sample(bones, size = 2, replace = TRUE)
    sum(dice)  # return last statement's return value
}

roll()

my_plot <- function(x, y, type = "l", ...) {
    plot(x, y, type = type, ...) # Pass '...' to 'plot' function
}
my_plot(x = 1L:5L, y = 2L:6L, main = "my plot")
