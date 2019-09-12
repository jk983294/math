library(plm)

time <- c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2)
firm <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6)
x <- c(1, 2, 4, 3, 5, 4, 2, 5, 5, 6, 7, 8)
y <- c(2, 3, 5, 4, 3, 2, 4, 7, 7, 8, 5, 7)
p.df <- data.frame(time, firm, x, y)
p.df <- pdata.frame(p.df, index = c("firm", "time"), drop.index = F, row.names = T)

# Clustering
pooled.plm <- plm(formula = y ~ x, data = p.df, model = "pooling")
coeftest(pooled.plm, vcov = vcovHC(pooled.plm, type = "sss", cluster = "group"))  # clustered by group
coeftest(pooled.plm, vcov = vcovHC(pooled.plm, type = "sss", cluster = "time"))  # clustered by time
coeftest(pooled.plm, vcov = vcovDC(pooled.plm, type = "sss"))  # Two-way clustering of both time and group

# Fixed Effect Model
fe.plm <- plm(formula = y ~ x, data = p.df, model = "within")
coeftest(fe.plm, vcov = vcovHC(fe.plm, type = "sss", cluster = "group"))
coeftest(fe.plm, vcov = vcovHC(fe.plm, type = "sss", cluster = "time"))
