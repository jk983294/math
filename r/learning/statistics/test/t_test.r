library(MASS)
# 针对两组的独立样本t检验可以用于检验两个总体的均值相等的假设
t.test(Prob ~ So, data = UScrime)  # p-value = 0.0006506, reject null hypothesis, mean is not equal
