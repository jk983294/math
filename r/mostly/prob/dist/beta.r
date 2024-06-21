# draw beta(a, b) density function
x <- (0L:100L) * 0.01
curve(dbeta(x, shape1 = 1., shape2 = 1.), from = 0., to = 1.)
curve(dbeta(x, shape1 = 2., shape2 = 1.), from = 0., to = 1.)
curve(dbeta(x, shape1 = 0.5, shape2 = 0.5), from = 0., to = 1.)
curve(dbeta(x, shape1 = 2, shape2 = 2), from = 0., to = 1.)

# often used as prior for a param in (0, 1)
# "conjugate prior to Binomial"
# X|p ~ Bin(n, p), p ~ Beta(a, b) [prior]
# find posterior dist p|X
# f(p|X=k) = p(X=k|p)*f(p) / P(X=k)
# p|X ~ Beta(a + X, b + n - X), X是观测到的成功次数
nsim <- 10^5
n <- 10
p <- runif(nsim) # p ~ Beta(1, 1) [prior]
x <- rbinom(nsim, n, p)
hist(x, breaks = seq(-0.5, n + 0.5, 1))
hist(p[x == 3]) # p|X=3 ~ Beta(4, 8) [posterior]

# 故事: X~Γ(a, λ), Y~Γ(b, λ)
# T= X+Y ~ Γ(a + b, λ) 总等待时间等于(a+b)次iid指数分布
# W=X/(X+Y) ~ Beta(a, b) 表示第一次等待时间占总等待时间的比例
# T和W独立
