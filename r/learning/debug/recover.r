f <- function(x, y) {
    z <- x + y
    g(z)
}
g <- function(x) {
    z <- round(x)
    h(z)
}

h <- function(x) {
    set.seed(1234)
    z <- rnorm(x)
    print(z)
}

options(error = recover)  # enable debugging mode

f(2, 3)
f(2, -3)  # enters debugging mode at this point

## type
## number - choose frame
## c - quit from that frame and choose another frame
## q - help

options(error = NULL)  # disable debugging mode
