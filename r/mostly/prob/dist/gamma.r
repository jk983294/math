# Γ(n+1) = n*Γ(n) = n!
# Γ(1/2) = sqrt(pi), Γ(1) = 1
# Γ(a, λ) 表示n个iid指数分布Expo(λ)的和, 总等待时间
# pdf: f(x) = (1/Γ(a)) * (λx)^(a-1) * e^(-λx) * (1/x)
x <- (0L:20L) * 0.01
curve(dgamma(x, shape = 3., rate = 1.), from = 0., to = 20.)
curve(dgamma(x, shape = 10., rate = 1), from = 0., to = 20.)
curve(dgamma(x, shape = 3., rate = 0.5), from = 0., to = 20.)
curve(dgamma(x, shape = 5., rate = 0.5), from = 0., to = 20.)
