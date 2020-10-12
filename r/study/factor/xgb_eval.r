library(data.table)
library(corrplot)
library(corrr)
library(dplyr)

dt1 <- fread("/tmp/eval.xgb.csv")
dim(dt1)

cor(dt1)
dt_mat <- as.matrix(dt1)
res <- rcorr(dt_mat)
res

corrplot(res$r, type = "upper", order = "hclust", p.mat = res$P)

draw_corr_bar_group_by <- function(dt, field) {
  cor_by_ <- dt %>% select(as.symbol(field), ret, ret_hat) %>% group_by(!!as.symbol(field)) %>% do({ correlate(select(., ret:ret_hat)) }) %>% filter(!is.na(ret)) %>% select(as.symbol(field), ret)
  ggplot(data=cor_by_, aes(x=!!as.symbol(field), y=ret)) + geom_bar(stat="identity")
}

# 按pi画出 ret/ret_hat的 corr
draw_corr_bar_group_by(dt1, "pi")

# 按di画出 ret/ret_hat的 corr
draw_corr_bar_group_by(dt1, "di")

# 按ti画出 ret/ret_hat的 corr
draw_corr_bar_group_by(dt1, "ti")
