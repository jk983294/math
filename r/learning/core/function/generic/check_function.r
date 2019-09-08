summary(women)
fit <- lm(weight ~ height, data = women)
summary(fit)

class(women)  # 'data.frame'
class(fit)  # 'lm'
methods(summary)

# see code of visible function
summary.data.frame
summary.lm
summary.default

# see code of 不可见的函数(在方法列表中加星号的函数)
getAnywhere(summary.ecdf)
