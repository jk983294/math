source("utils/utils.r")

dat <- read_stock_data()
dr <- dailyReturn(dat)
my_stats <- mystats(dr$daily.returns)
kurtosis <- my_stats["kurtosis"]

# t test
t1 <- kurtosis/sqrt(24/my_stats["n"])
pvalue <- 2 * (1 - pnorm(t1))
if (pvalue < 0.05) {
    print("kurtosis excess != 0, fat tail")
} else {
    print("kurtosis excess == 0")
}
