# on vector
price <- 1:10
initial <- 5
x <- Reduce(function(x, y) {
    (x + 2 * y)/3
}, price, initial)
x

# on list
HI <- data.frame(stock = c(1, 2), high = c(7, 8))
LO <- data.frame(stock = c(1, 2), low = c(1, 2))
OP <- data.frame(stock = c(1, 2), open = c(3, 6))
CL <- data.frame(stock = c(1, 2), close = c(5, 5))
df <- Reduce(function(x, y) {
    merge(x, y)  # since merging is one by one on the data.frame, we just merge all in list
}, list(HI, LO, OP, CL))
df
