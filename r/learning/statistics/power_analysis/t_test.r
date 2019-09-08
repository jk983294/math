library(pwr)
pwr.t.test(d = 0.8, sig.level = 0.05, power = 0.9, type = "two.sample", alternative = "two.sided")
pwr.t.test(n = 20, d = 0.5, sig.level = 0.01, type = "two.sample", alternative = "two.sided")
