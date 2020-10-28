library(foreach)
library(doParallel)
registerDoParallel(8)

# 顺序版本
foreach(i = 1:3) %do% {
    sqrt(i)
}

# 并行版本
foreach(i = 1:3) %dopar% {
    sqrt(i)
}

# 并行版本，Return a vector
foreach(i = 1:3, .combine = c) %dopar% {
    sqrt(i)
}

# 并行版本，closure to set result
data <- 4:10
result <- foreach(i = seq(data)) %dopar% {
    my_ret_list <- list(d = sqrt(data[i]))
    return(my_ret_list)
}

vec_result <- sapply(result, function(x) x$d)
