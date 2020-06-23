dyn.load("./library/rcpp/compile/fibonacci.so")
.Call("fibWrapper", 5)
