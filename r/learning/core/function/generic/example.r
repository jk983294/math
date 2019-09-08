mymethod <- function(x, ...) UseMethod("mymethod")
mymethod.a <- function(x) print("Using A")
mymethod.b <- function(x) print("Using B")
mymethod.default <- function(x) print("Using Default")

x <- 1:5
y <- 6:10
z <- 10:15
class(x) <- "a"
class(y) <- "b"

mymethod(x) # "Using A"
mymethod(y) # "Using B"
mymethod(z) # "Using Default"

class(z) <- c("a", "b")
mymethod(z) # "Using A"  第一类用来决定哪个泛型函数被调用
class(z) <- c("c", "a", "b")
mymethod(z) # "Using A"  没有mymethod.c()函数,因此下一个类"a"被使用
