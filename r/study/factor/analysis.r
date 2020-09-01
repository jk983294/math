library(data.table)
library(Hmisc)
library(corrplot)
library(PerformanceAnalytics)

dt1 <- fread("~/github/MyTmp/sim/cmake-build-debug/pmfut/factor.csv")
dim(dt1)

cor(dt1)
dt_mat <- as.matrix(dt1)
res <- rcorr(dt_mat)
res

corrplot(res$r, type = "upper", order = "hclust", p.mat = res$P)
# chart.Correlation(dt_mat, histogram=FALSE, pch=19)
