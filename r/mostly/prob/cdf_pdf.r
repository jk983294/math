# 求连续随机变量quantile
# CDF的根等于quantile
q_ <- 0.5
g <- function(x) pexp(x) - q_
uniroot(g, lower = 0., upper = 1.)

# 求连续随机变量mode
# 就是求pdf的最大值对应的x
h <- function(x) x^5 * exp(-x)
optimize(h, lower = 0, upper = 20, maximum = TRUE)

# 求离散随机变量quantile
# CDF的根大于等于quantile
n <- 50
which.max(pbinom(0:n, n, 0.2) >= q_)

# 求离散随机变量mode
# 就是求pdf的最大值对应的x
which.max(dbinom(0:n, n, 0.2))
