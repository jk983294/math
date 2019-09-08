trim <- function(p) {
    trimit <- function(x) {
        n <- length(x)
        lo <- floor(n * p) + 1
        hi <- n + 1 - lo
        x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
    }
    trimit
}
x <- 1:10
trim10pct <- trim(0.1)
y <- trim10pct(x)
y  # [2 3 4 5 6 7 8 9]
trim20pct <- trim(0.2)
y <- trim20pct(x)
y  # [3 4 5 6 7 8]

# argument bind to function env
ls(environment(trim10pct))
get("p", env = environment(trim10pct))  # 0.1
ls(environment(trim20pct))
get("p", env = environment(trim20pct))  # 0.2

# another example
makeFunction <- function(k) {
    f <- function(x) {
        print(x + k)
    }
}

g <- makeFunction(10)
g(4)  # 14
k <- 2
g(5)  # 15

ls(environment(g))
environment(g)$k  # 10
