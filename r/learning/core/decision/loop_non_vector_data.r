# loop non-vector data
u <- 1
v <- 2
for (m in c("u", "v")) {
    z <- get(m)
    print(z)
}
