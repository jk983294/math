src <- "
    Function sort(x) ;
    return sort( y, Named(\"decreasing\", true));
"

fun <- cxxfunction(signature(x = "function", y = "ANY"), src, plugin = "Rcpp")

# pass in functor sort
fun(sort, sample(1:5, 10, TRUE))
fun(sort, sample(LETTERS[1:5], 10, TRUE))


src <- "
    RNGScope scp;
    Rcpp::Function rt(\"rt\"); // access R function by string name
    return rt(5, 3);
"

fun <- cxxfunction(signature(), src, plugin = "Rcpp")
set.seed(42)
fun()
fun()

set.seed(42)
rt(10, 3)  # draw 10 number from t-distribution with dof 3
