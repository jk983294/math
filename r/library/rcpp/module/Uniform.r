setClass("Uniform", representation(pointer = "externalptr"))

# helper
Uniform_method <- function(name) {
    paste("Uniform", name, sep = "__")
}

# syntactic sugar to allow object$method( ... )
setMethod("$", "Uniform", function(x, name) {
    function(...) .Call(Uniform_method(name), x@pointer, ...)
})

# syntactic sugar to allow new( 'Uniform', ... )
setMethod("initialize", "Uniform", function(.Object, ...) {
    .Object@pointer <- .Call(Uniform_method("new"), ...)
    .Object
})

u <- new("Uniform", 0, 10)
u$draw(10L)
