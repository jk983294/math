# frame是一种特殊的列表,集合中每个原子向量都有相同的长度
head(iris)
unclass(iris)  # 每个向量代表数据框中的一列
attributes(iris)  # names属性(column name), row.names属性和class属性'data.frame'
