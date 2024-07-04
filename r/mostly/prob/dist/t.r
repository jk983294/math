# useful in situations where the sample size is small and the population variance is unknown
# it is similar in shape to the standard normal distribution but has heavier tails
# T = Z / sqrt(V / n) where Z ~ N(0, 1), V ~ χ^2(n)
# t-dist converge to std normal dist when n go large
x <- (-10L:10L) * 0.1
curve(dt(x, df = 5), from = -10., to = 10)
curve(dt(x, df = 10), from = -10., to = 10)
curve(dt(x, df = 20), from = -10., to = 10)
