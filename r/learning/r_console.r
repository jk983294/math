# install package
install.packages("ggplot2", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
install.packages("formatR", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
install.packages("reshape", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
install.packages("plotrix", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

# help
help.start()  # 打开帮助文档首页
help("foo")  # ?foo 查看函数 foo 的帮助(引号可以省略)
help(package = "ggplot2")  # package desc
help.search("foo")  # ??foo 以 foo 为关键词搜索本地帮助文档
example("foo")  # 函数 foo 的使用示例
RSiteSearch("foo")  # 以 foo 为关键词搜索在线文档和邮件列表存档
apropos("foo", mode = "function")  # 列出名称中含有 foo 的所有可用函数
data()  # 列出当前已加载包中所含的所有可用示例数据集
vignette()  # 列出当前已安装包中所有可用的 vignette 文档
vignette("foo")  # 为主题 foo 显示指定的 vignette 文档
args(round)  # check args of a given function
