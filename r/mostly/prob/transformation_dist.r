# 单变量 Y = g(X)
# f_Y(y) = f_X(x) * |dx/dy|
# 多元变量 Y = g(X)
# f_Y(y) = f_X(x) * |det(dx/dy)| 其中dx/dy是Jacobian matrix, f_X(x)是X的联合分布, f_Y(y)是Y的联合分布

# Sum and Convolution
# T = X + Y, X~U(0,1), Y~U(0,1), X and Y iid
x <- runif(10^5)
y <- runif(10^5)
hist(x + y)
