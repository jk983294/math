## symbol and line
* pch - 指定绘制点时使用的符号
* cex - 指定符号的大小, cex是一个数值,表示绘图符号相对于默认大小的缩放倍数
* lty - 指定线条类型(参见图 3-5)
* lwd - 指定线条宽度, lwd是以默认值的相对大小来表示的(默认值为1)

Ex: plot(dose, drugA, type="b", lty=3, lwd=3, pch=15, cex=2)

## color
* col - 默认的绘图颜色。某些函数(如 lines和pie)可以接受一个含有颜色值的向量并自动循环使用
* col.axis - 坐标轴刻度文字的颜色
* col.lab - 坐标轴标签(名称)的颜色
* col.main - 标题颜色
* col.sub - 副标题颜色
* fg - 图形的前景色
* bg - 图形的背景色

white: col=1, col="white", col="#FFFFFF", col=rgb(1,1,1), col=hsv(0,0,1)

## text
* cex - 表示相对于默认大小缩放倍数的数值
* cex.axis - 坐标轴刻度文字的缩放倍数
* cex.lab - 坐标轴标签(名称)的缩放倍数
* cex.main - 标题的缩放倍数
* cex.sub - 副标题的缩放倍数
* font - 整数, 用于指定绘图使用的字体样式, 1=常规, 2=粗体, 3=斜体, 4=粗斜体, 5=符号字体
* font.axis - 坐标轴刻度文字的字体样式
* font.lab - 坐标轴标签(名称)的字体样式
* font.main - 标题的字体样式
* font.sub - 副标题的字体样式
* ps - 字体磅值, 文本的最终大小为 ps * cex
* family - 绘制文本时使用的字体族, 标准的取值为serif(衬线), sans(无衬线) 和mono(等宽)

Ex: par(font.lab=3, cex.lab=1.5, font.main=4, cex.main=2)

## 图形尺寸与边界尺寸
* pin - 以英寸表示的图形尺寸(宽和高)
* mai - 以数值向量表示的边界大小, 顺序为“下、左、上、右”, 单位为英寸
* mar - 以数值向量表示的边界大小, 顺序为“下、左、上、右”, 单位为英分

## title
```R
title(main="My Title", col.main="red",
        sub="My Subtitle", col.sub="blue",
        xlab="My X label", ylab="My Y label",
        col.lab="green", cex.lab=0.75)
```
