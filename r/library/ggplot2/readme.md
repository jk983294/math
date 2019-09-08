## overview
aes()函数的功能是指定每个变量扮演的角色 (aes代表aesthetics,即如何用视觉形式呈现信息)

使用一个或多个几何函数向图中添加了几何对象(简写为geom),包括点、线、条、箱线图和阴影区域

## type
| 函数        | 添加   |  选项  |
| --------   | -----:  | :----:  |
| geom_bar()        |  条形图  |   color 、 fill 、 alpha  |
| geom_boxplot()    |  箱线图  |   color 、 fill 、 alpha 、 notch 、 width  |
| geom_density()    |  密度图  |   color 、 fill 、 alpha 、 linetype  |
| geom_histogram()  |  直方图  |   color 、 fill 、 alpha 、 linetype 、 binwidth  |
| geom_hline()      |  水平线  |   color 、 alpha 、 linetype 、 size  |
| geom_jitter()     |  抖动点  |   color 、 size 、 alpha 、 shape  |
| geom_line()       |  线图   |    colorvalpha 、 linetype 、 size  |
| geom_point()      |  散点图  |   color 、 alpha 、 shape 、 size  |
| geom_rug()        |  地毯图  |   color 、 side  |
| geom_smooth()     |  拟合曲线 |   method 、 formula 、 color 、 fill 、 linetype 、 size  |
| geom_text()       |  文字注解 |  |
| geom_violin()     |  小提琴图 |   color 、 fill 、 alpha 、 linetype  |
| geom_vline()      |  垂线   |    color 、 alpha 、 linetype 、 size  |

## options
* color     -  对点、线和填充区域的边界进行着色
* fill      -  对填充区域着色,如条形和密度区域
* alpha     -  颜色的透明度,从0(完全透明)到1(不透明)
* linetype  -  图案的线条(1=实线,2=虚线,3=点,4=点破折号,5=长破折号,6=双破折号)
* size      -  点的尺寸和线的宽度
* shape     -  点的形状(和pch一样,0=开放的方形,1=开放的圆形,2=开放的三角形,等等)
* position  -  绘制诸如条形图和点等对象的位置。对条形图来说, "dodge" 将分组条形图并排, "stacked" 堆叠分组条形图, "fill" 垂直地堆叠分组条形图并规范其高度相等。对于点来说, "jitter" 减少点重叠
* binwidth  -  直方图的宽度
* notch     -  表示方块图是否应为缺口(TRUE/FALSE)
* sides     -  地毯图的安置( "b" =底部, "l" =左部, "t" =顶部, "r" =右部, "bl" =左下部,等等)
* width     -  箱线图的宽度

## smooth
* method     - 使用的平滑函数。允许的值包括 lm 、 glm 、 smooth 、 rlm 和 gam ,分别对应线性、广义线性、 loess 、健壮线性和广义相加模型。 smooth 是默认值
* formula    - 在光滑函数中使用的公式。例子包括y~x(默认), y~log(x),y~poly(x,n) 表示n次多项式拟合y~ns(x,n) 表示一个具有n个自由度的样条拟合
* se         - 绘制置信区间(TRUE/FALSE)。默认为TRUE
* level      - 使用的置信区间水平(默认为95%)
* fullrange  - 指定拟合应涵盖全图(TRUE)或仅仅是数据(FALSE)。 默认为FALSE

## axis
