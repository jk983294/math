## PCA overview
主成分分析(PCA)是一种数据降维技巧,它能将大量相关变量转化为一组很少的不相关变量,这些无关变量称为主成分

主成分(PC1和PC2)是观测变量(X1到X5)的线性组合。形成线性组合的权重都是通过最大化各主成分所解释的方差来获得,同时还要保证各主成分间不相关

* PC1 栏包含了成分载荷,指观测变量与主成分的相关系数。成分载荷(component loadings)可用来解释主成分的含义。
* h2 栏指成分公因子方差,即主成分对每个变量的方差解释度。
* u2 栏指成分唯一性,即方差无法被主成分解释的比例(1–h2)。
* SS loadings 行包含了与主成分相关联的特征值,指的是与特定主成分相关联的标准化后的方差值
* Proportion Var 行表示的是每个主成分对整个数据集的解释程度。
* Cumulative Var - Cumulative of Proportion Var

## PCA rotate
旋转是一系列将成分载荷阵变得更容易解释的数学方法,它们尽可能地对成分去噪。

旋转方法有两种:使选择的成分保持不相关(正交旋转),和让它们变得相关(斜交旋转)。

旋转方法也会依据去噪定义的不同而不同。

最流行的正交旋转是方差极大旋转,它试图对载荷阵的列进行去噪,使得每个成分只由一组有限的变量来解释(即载荷阵每列只有少数几个很大的载荷,其他都是很小的载荷)。

## EFA overview
EFA是一系列用来发现一组变量的潜在结构的方法。它通过寻找一组更小的、潜在的或隐藏的结构来解释已观测到的、显式的变量间的关系

这些虚拟的、无法观测的变量称作因子。

对于EFA,Kaiser-Harris准则的特征值数大于0,而不是1。

## power analysis
两种方法都需要大样本来支撑稳定的结果

因子分析需要5~10倍于变量数的样本数
