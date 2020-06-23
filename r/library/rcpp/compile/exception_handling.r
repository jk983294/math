dyn.load("./library/rcpp/compile/fibonacci.so")
.Call("fun", 4)
.Call("fun", -4)
.Call("fun", 11)
