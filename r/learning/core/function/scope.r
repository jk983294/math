# 在函数之外创建的对象是全局的(也适用于函数内部)。在函数之内创建的对象是局部的(仅仅适用于函数内部).
# 局部对象在函数执行后被丢弃。只有那些通过return()(或使用<<-分配)传回的对象在函数执行之后可以继续使用.
# 全局对象在函数之内可被访问(可读)但是不会改变(除非使用<<-).
# 对象可以通过参数传递到函数中,但是不会被函数改变。传递的是对象的副本而不是变量本身.

x <- 2
y <- 3
z <- 4
f <- function(w) {
    x <- w * y * z
    z <<- 2
    assign("y", 1, pos = .GlobalEnv)
    return(x)
}
f(x)  # 12
x
y  # changed to 1 since we use assign inside function
z  # changed to 2 since we use <<- inside function

# Lexical Scoping
y <- 10L
f <- function(x) {
    y <- 2L
    y^2L + g(x)
}
g <- function(x) {
    x * y
}
f(3L) # y in g is 10, not 2 defined in f
