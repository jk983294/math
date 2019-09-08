# independence in contingency tables
library(coin)
library(vcd)
Arthritis <- transform(Arthritis, Improved = as.factor(as.numeric(Improved)))
set.seed(1234)
# check if Treatment and Improved are independent
chisq_test(Treatment ~ Improved, data = Arthritis, distribution = approximate(nresample = 9999))

# independence between numeric variables
states <- as.data.frame(state.x77)
set.seed(1234)
# check if Illiteracy and Murder are independent
spearman_test(Illiteracy ~ Murder, data = states, distribution = approximate(nresample = 9999))
