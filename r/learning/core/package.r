# install package
install.packages("ggplot2", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
install.packages("formatR", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
install.packages("vcd", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")  # Visualizing Categorical Data
install.packages("RColorBrewer")  # color set
install.packages("vioplot")  # 小提琴图
install.packages("gplots")  # 带有置信区间的组均值图形
install.packages("scatterplot3d")
install.packages("coefplot")  # coefficient plot for the VAR model
# sudo apt-get install xorg libx11-dev libglu1-mesa-dev # rgl dependency
install.packages("rgl")  # 3d plot
install.packages("corrgram")  # correlogram
install.packages("ggthemes")
install.packages("dygraphs")

install.packages("Hmisc")
install.packages("reshape2")
install.packages("sm")  # smooth method
# sudo apt-get install libgmp-dev libmpfr-dev # HH dependency
install.packages("HH")  # 绘制因变量、协变量和因子之间的关系图
install.packages("pastecs")  # statistics
install.packages("psych")  # statistics
install.packages("doBy")  # statistics
install.packages("gmodels")  # cross table
# sudo apt-get install libblas-dev liblapack-dev # ggm dependency
install.packages("ggm")  # 偏相关系数
# sudo apt-get install libcurl4-openssl-dev # car dependency
install.packages("car", dependencies = TRUE)  # scatterplot
install.packages("effects")  # 图形展示交互项的结果
install.packages("leaps")  # 全子集回归
install.packages("bootstrap")  # k重交叉验证
install.packages("multcomp")  # 多重均值
install.packages("pwr")  # Power analysis
install.packages("coin")  # 置换检验的框架
install.packages("lmPerm")  # 方差分析和回归分析的置换检验
install.packages("qcc")  # 泊松模型过度离势的检验方法
install.packages("psych")  # 主成分分析, 因子分析
install.packages("GPArotation")  # 因子分析
install.packages("forecast")  # 时序分析
install.packages("rugarch")  # ARCH and GARCH model
install.packages("vars")  # VAR model
install.packages(c("VIM", "mice"))  # missing data
install.packages("xts")  # time series
install.packages(c("sandwich", "lmtest"))  # robust standard error for lm

# core
install.packages("doParallel", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
install.packages("roxygen2", depend = TRUE)
# sudo apt install texlive-latex-base texlive-latex-recommended texlive-fonts-recommended texlive-latex-extra


# data
install.packages("AER")
install.packages("robust")

# install from local file
install.packages("~/XML_3.98-1.3.zip", repos = NULL, type = "source")

# import package
library("ggplot2")
library("package Name", lib.loc = "path to library")

# update packages
update.packages()

# location /usr/local/lib/R/site-library
.libPaths()  # get library locations
library()  # get all packages installed
search()  # all packages currently loaded
installed.packages()  # list installed packages metadata
