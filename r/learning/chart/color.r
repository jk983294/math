colors()  # 返回所有可用颜色的名称

# create color vector
rainbow(5)
heat.colors(5)
terrain.colors(5)
topo.colors(5)
cm.colors(5)

# create color vector from RColorBrewer
library(RColorBrewer)
n <- 7
mycolors <- brewer.pal(n, "Set1")
barplot(rep(1, n), col = mycolors)
brewer.pal.info  # 得到所有可选调色板的列表, like Set1

# 多阶灰度色
n <- 10
mycolors <- rainbow(n)
pie(rep(1, n), labels = mycolors, col = mycolors)
mygrays <- gray(0:n/n)
pie(rep(1, n), labels = mygrays, col = mygrays)
