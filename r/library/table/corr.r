library(data.table)
library(corrr)
library(dplyr)

ii <- c(1,2,1,2)
Time <- c(1,1,2,2)
signal1 <- c(1.1,2.2,3.3,4.4)
signal2 <- c(1.2,2.3,3.4,4.5)
signals <- data.table(Time,ii,signal1, signal2)

# signal1, signal2 整列corr
signals %>% select(signal1, signal2) %>% correlate()

# 按II做signal1, signal2列的corr
signals %>% group_by(ii) %>% do({ correlate(select(., signal1:signal2)) }) # cor group by ii

# draw
x <- signals %>% select(signal1, signal2) %>% correlate() %>% rearrange() %>% shave()
fashion(x)
rplot(x)
