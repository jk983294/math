# It's' the distribution of the ratio of two independent normally distributed random variables with mean zero.
# X/Y ~ Cauchy(0, 1)
# its expected value and its variance are undefined
# The Cauchy distribution has no moment generating function
# it has fat tail against normal dist
hist(rcauchy(10^3))
x <- (-100L:100L) * 0.1
curve(dcauchy(x, location = 0., scale = 1.), from = -10., to = 10.)
