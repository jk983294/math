# Chi-square test of independence 对二维表的行变量和列变量进行卡方独立性检验
library(vcd)
mytable <- xtabs(~Treatment + Improved, data = Arthritis)
chisq.test(mytable)  # p-value = 0.001463 less than the significance level (0.05), then reject null hypothesis, so they are not independent
mytable <- xtabs(~Improved + Sex, data = Arthritis)
chisq.test(mytable)  # p-value = 0.08889 > 0.05, they are independent

# Fisher's exact test
mytable <- xtabs(~Treatment + Improved, data = Arthritis)
fisher.test(mytable)  # p-value = 0.001393, they are not independent

# Chochran-Mantel-Haenszel test
mytable <- xtabs(~Treatment + Improved + Sex, data = Arthritis)
mantelhaen.test(mytable)
