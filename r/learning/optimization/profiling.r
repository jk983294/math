Rprof("profile1.out", line.profiling = TRUE)

# create big data frame:
n <- 1000
x <- data.frame(group = sample(letters[1:4], n, replace = TRUE), condition = sample(LETTERS[1:10], 
    n, replace = TRUE), data = rnorm(n))

# reasonable operations:
marginal.means.1 <- aggregate(data ~ group + condition, data = x, FUN = mean)

# unreasonable operations:
marginal.means.2 <- marginal.means.1[NULL, ]

row.counter <- 1
for (condition in levels(x$condition)) {
    for (group in levels(x$group)) {
        tmp.value <- 0
        tmp.length <- 0
        for (c in 1:nrow(x)) {
            if ((x[c, "group"] == group) & (x[c, "condition"] == condition)) {
                tmp.value <- tmp.value + x[c, "data"]
                tmp.length <- tmp.length + 1
            }
        }
        marginal.means.2[row.counter, "group"] <- group
        marginal.means.2[row.counter, "condition"] <- condition
        marginal.means.2[row.counter, "data"] <- tmp.value/tmp.length
        row.counter <- row.counter + 1
    }
}

# does it produce the same results?
all.equal(marginal.means.1, marginal.means.2)

Rprof(NULL)

summaryRprof("profile1.out", lines = "show")

library(profr)
plot(parse_rprof("profile1.out"))
