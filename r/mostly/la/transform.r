# point
point <- matrix(c(1, 1), nrow = 2)
points <- matrix(c(1, 1, 1, 3), nrow = 2)

# transform matrix
m0 <- matrix(c(2, 0, 0, 3), nrow = 2) # scale
m0 %*% points # 按行scale
points %*% m0 # 按列scale

m1 <- matrix(c(1, 0, 0, -1), nrow = 2) # 关于x轴的对称
m2 <- matrix(c(-1, 0, 0, 1), nrow = 2) # 关于y轴的对称
m3 <- matrix(c(0, 1, 1, 0), nrow = 2) # 关于y=x的对称
m4 <- matrix(c(0, -1, -1, 0), nrow = 2) # 关于y=-x的对称
m5 <- matrix(c(-1, 0, 0, -1), nrow = 2) # 关于(0, 0)的对称
m6 <- matrix(c(0.5, 0, 0, 1), nrow = 2) # 水平收缩拉伸
m7 <- matrix(c(1, 0, 0, 1.5), nrow = 2) # 垂直收缩拉伸
m8 <- matrix(c(1, 0, 0.5, 1), nrow = 2) # 水平剪切
m9 <- matrix(c(1, 1.5, 0, 1), nrow = 2) # 垂直剪切
m10 <- matrix(c(1, 0, 0, 0), nrow = 2) # 水平投影
m11 <- matrix(c(0, 0, 0, 1), nrow = 2) # 垂直投影
angle <- pi / 4.
m12 <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), nrow = 2) # 逆时针旋转pi/4
m13 <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow = 2) # 顺时针旋转pi/4
move_func <- function(point, delta_x, delta_y) { # 二维平面平移
    move_mat <- matrix(c(1, 0, 0, 0, 1, 0, delta_x, delta_y, 1.), 3, 3)
    point_ <- c(point, 1.) # 齐次坐标
    moved <- move_mat %*% point_
    moved[0:2]
}

m1 %*% point
m2 %*% point
m3 %*% point
m4 %*% point
m5 %*% point
m6 %*% point
m7 %*% point
m8 %*% point
m9 %*% point
m10 %*% point
m11 %*% point
m12 %*% point
m13 %*% point
move_func(point, 1., 2.)
