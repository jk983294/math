# testing a correlation coefficient for significance.
# 当研究的假设为总体的相关系数小于0时, 使用alternative='less'.
# 在研究的假设为总体的相关系数大于0时, 使用alternative='greater'.
# 在默认情况下,假设为alternative='two.side' (相关系数不等于0).

# The P value is the probability of finding the observed results when the null
# hypothesis (H0) is true. the lower the p-value, the more ridiculous our null
# hypothesis looks. If the p-value is lower than a predetermined significance
# level, then we reject the null hypothesis.
cor.test(states[, 3], states[, 5])


# correlation matrix and tests of significance via corr.test
library(psych)
corr.test(states, use = "complete")
