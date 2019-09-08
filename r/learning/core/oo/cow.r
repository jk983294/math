a <- 1
z <- a

tracemem(a)
tracemem(z) # a z same address
a <- 2
tracemem(a) # a address changed
tracemem(z) # z address remain

z <- runif(5)
tracemem(z)
z[3] <- 8
tracemem(z) # vector address won't change if element re-assigned
