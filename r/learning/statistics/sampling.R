############# sample default is without replacement
sample(x = 1:4, size = 2)  # pick size element from x

# roll a dice
die <- 1:6
sample(x = die, size = 2)
sample(x = die, size = 2, replace = TRUE)  # with replacement sampling

# experiment many times
roll <- function(bones = 1:6) {
    dice <- sample(bones, size = 2, replace = TRUE)
    sum(dice)  # return last statement's return value
}
replicate(10, roll())
