# uniformly distributed random numbers, [0, 1] by default
runif(3)
runif(3, min = 0, max = 100)
floor(runif(3, min = 0, max = 101))  # [0, 100] integer
floor(runif(3, min = 0, max = 2))  # [0, 1] integer

# normal distribution
rnorm(3)
rnorm(3, mean = 2, sd = 5)

# seed 让结果可以重现
set.seed(1234)
