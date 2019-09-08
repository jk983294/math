cor(mtcars)

library(corrgram)

# 蓝色和从左下指向右上的斜杠表示单元格中的两个变量呈正相关。反过来,红色和从左上指向右下的斜杠表示变量呈负相关。
# 色彩越深,饱和度越高,说明变量相关性越大
corrgram(mtcars, order = TRUE, lower.panel = panel.shade, upper.panel = panel.pie, 
    text.panel = panel.txt, main = "Corrgram of mtcars intercorrelations")

corrgram(mtcars, order = TRUE, lower.panel = panel.ellipse, upper.panel = panel.pts, 
    text.panel = panel.txt, diag.panel = panel.minmax, main = "Corrgram of mtcars data using scatter plots
and ellipses")

# use other color system
cols <- colorRampPalette(c("darkgoldenrod4", "burlywood1", "darkkhaki", "darkgreen"))
corrgram(mtcars, order = TRUE, col.regions = cols, lower.panel = panel.shade, upper.panel = panel.conf, 
    text.panel = panel.txt, main = "A Corrgram (or Horse) of a Different Color")
