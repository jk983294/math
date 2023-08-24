(x1 <- 0:10)
(x2 <- 10:0) # Descending order
(x3 <- seq(10))
(x4 <- seq(30, 0, by = -3))
(x5 <- c(5, 4, 1, 6, 7, 2, 2, 3, 2, 8))
(x6 <- rep(TRUE, 5))
(x7 <- rep(c(TRUE, FALSE), 5))
(x8 <- rep(c(TRUE, FALSE), each = 5))
(x8 <- rep(letters[1:3], each = 2))
(x8 <- rep(LETTERS[1:3], each = 2))

(x1 <- runif(5)) # [0, 1]
(x1 <- runif(5, min = 0, max = 100)) # [0, 100]
(x1 <- floor(runif(5, min = 0, max = 100))) # random int
(x1 <- rnorm(5))
(x1 <- rnorm(5, mean=5, sd=10))