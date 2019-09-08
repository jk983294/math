library("ggplot2")

x <- c(0, 0, 0, 0, 1, 2, 2, 2, 3, 3)
qplot(x, binwidth = 1)

roll <- function() {
    # roll a unfair dice
    die <- 1:6
    prob <- c(1/8, 1/8, 1/8, 1/8, 1/8, 3/8)
    dice <- sample(die, size = 2, replace = TRUE, prob = prob)
    sum(dice)  # return last statement's return value
}
rolls <- replicate(1e+05, roll())
qplot(rolls, binwidth = 1)
