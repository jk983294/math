library(Rcpp)

fibR <- function(n) {
    if (n == 0) 
        return(0)
    if (n == 1) 
        return(1)
    return(fibR(n - 1) + fibR(n - 2))
}

## memoization solution courtesy of Pat Burns
mfibR <- local({
    memo <- c(1, 1, rep(NA, 1000))
    f <- function(x) {
        if (x == 0) 
            return(0)
        if (x < 0) 
            return(NA)
        if (x > length(memo)) 
            stop("x too big for implementation")
        if (!is.na(memo[x])) 
            return(memo[x])
        ans <- f(x - 2) + f(x - 1)
        memo[x] <<- ans
        ans
    }
})

## linear / iterative solution
fibRiter <- function(n) {
    first <- 0
    second <- 1
    third <- 0
    for (i in seq_len(n)) {
        third <- first + second
        first <- second
        second <- third
    }
    return(first)
}

fibR(5)
mfibR(5)
fibRiter(5)

sourceCpp("./library/rcpp/attribute/fibonacci.cpp")
fibonacci(5)

cpptxt <- "
int fibonacci(const int x) {
if (x < 2) return(x);
return (fibonacci(x - 1)) + fibonacci(x - 2);
}"
fibCpp <- cppFunction(cpptxt)
fibCpp(5)
