library(magrittr)

u <- function(x) x + 2.
f <- function(x) x * 3.
g <- function(x) x / 2.
h <- function(x, y) x + y - 1.

x <- c(1.:10.)

# piping
x %>% f # f(x)
x %>% h(1.) # f(x, 1.)
x %>% f %>% g %>% u # u(g(f(x)))

# argument placeholder .
x %>% h(1., .) # x %>% f(y, .) is equivalent to f(y, x)
x %>% h(1., y = .) # x %>% f(y, z = .) is equivalent to f(y, z = x)
# x %>% f(y = nrow(.), z = ncol(.)) is equivalent to f(x, y = nrow(x), z = ncol(x))
# x %>% {f(y = nrow(.), z = ncol(.))} is equivalent to f(y = nrow(x), z = ncol(x))

# build unary functions
f1 <- . %>% cos %>% sin # = f <- function(.) sin(cos(.))
data %>% f1

# %$%
# functions that do not have a data argument, for which it is useful to expose the variables in the data
iris %$% cor(Sepal.Length, Sepal.Width)
df <- data.frame(z = rnorm(100L))
df %$% ts.plot(z)

# assignment operator %<>%
df %<>% transform(z = z * 2.) # equals with mtcars <- mtcars %>% transform(cyl = cyl * 2)
