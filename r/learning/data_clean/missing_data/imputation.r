library(mice)
data(sleep, package = "VIM")
imp <- mice(sleep, seed = 1234) # 包含m个插补数据集的list
fit <- with(imp, lm(Dream ~ Span + Gest)) # 包含m个单独统计分析结果的list
pooled <- pool(fit) # m个统计分析平均结果的列表对象
summary(pooled)
imp
dataset3 <- complete(imp, action=3) # 观察插补数据集3
