t1 <- Sys.time()

for (i in 1:1e+05) {
    x <- x + i
}

t2 <- Sys.time()
t2 - t1
