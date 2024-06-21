# f 是一个凸函数 (f''(x) > 0)，X 是一个随机变量，then
# E[f(X)] ≥ f(E[X])
# 凸函数使得随机变量的平均值比随机变量函数值的平均值要小
# E(X^2) ≥ (EX)^2
# E(1/X) ≥ 1/EX
# E(lnX) ≤ ln(EX)

x <- rexp(10^4)
mean(log(x))
log(mean(x))
