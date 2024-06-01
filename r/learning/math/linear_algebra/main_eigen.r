A <- matrix(c(1.8, 0.2, 0.8, 1.2), nrow = 2)

get_main_eigen <- function(A, eps = 1e-6) {
    u_k <- 0.
    x_k <- rep(0., dim(A)[1])
    x_k[1] <- 2.
    loop_cnt <- 0
    repeat {
        loop_cnt <- loop_cnt + 1
        x_k <- A %*% x_k
        u_k_new <- max(abs(x_k))
        x_k <- x_k * (1. / u_k_new)
        if (abs(u_k_new - u_k) < eps) {
            break
        } else {
            u_k <- u_k_new
        }
    }
    list(eigen_val = u_k, eigen_vec = x_k, cnt = loop_cnt)
}
(res <- get_main_eigen(A))
