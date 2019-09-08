set.seed(1234)  # make results reproducible

# determine how many factors
library(psych)
# Kaiser-Harris准则建议保留特征值大于1的主成分,特征值小于1的成分所解释的方差比包含在单个变量中的方差更少
fa.parallel(USJudgeRatings[, -1], fa = "pc", n.iter = 100, show.legend = FALSE, main = "Scree plot with parallel analysis")  # 1 factor

# Principal components analysis of US Judge Ratings
library(psych)
pc <- principal(USJudgeRatings[, -1], nfactors = 1)
pc

# Principal components analysis Harman23.cor data
library(psych)
fa.parallel(Harman23.cor$cov, n.obs = 302, fa = "pc", n.iter = 100, show.legend = FALSE, 
    main = "Scree plot with parallel analysis")  # 2 factors

PC <- principal(Harman23.cor$cov, nfactors = 2, rotate = "none")
PC

# varimax rotation
rc <- principal(Harman23.cor$cov, nfactors = 2, rotate = "varimax")
rc


# Obtaining componenet scores, 每个观测在成分上的得分
pc <- principal(USJudgeRatings[, -1], nfactors = 1, score = TRUE)
head(pc$scores)
cor(USJudgeRatings$CONT, pc$score)  # 律师与法官的熟稔度与律师的评分毫无关联


# Obtaining principal component scoring coefficients
library(psych)
rc <- principal(Harman23.cor$cov, nfactors = 2, rotate = "varimax")
round(unclass(rc$weights), 2)
