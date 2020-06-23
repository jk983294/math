library(Rcpp)

cppFunction("
    int fun2(int dx) {
        if ( dx > 10 )
            throw std::range_error(\"too big\");
        return dx * dx;
    }
")

fun2(3)
fun2(13)
