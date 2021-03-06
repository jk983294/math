## overview
回归方法本就是用来从一堆数据中获取最优模型参数。
对于OLS回归,通过使得预测误差(残差)平方和最小和对响应变量的解释度(R平方)最大,可获得模型参数

## method
* 简单线性    - 用一个量化的解释变量预测一个量化的响应变量
* 多项式      - 用一个量化的解释变量预测一个量化的响应变量,模型的关系是 n 阶多项式
* 多层        - 用拥有等级结构的数据预测一个响应变量(例如学校中教室里的学生)。也被称为分层模型、嵌套模型或混合模型
* 多元线性     - 用两个或多个量化的解释变量预测一个量化的响应变量
* 多变量       - 用一个或多个解释变量预测多个响应变量
* Logistic    - 用一个或多个解释变量预测一个类别型响应变量
* 泊松        - 用一个或多个解释变量预测一个代表频数的响应变量
* Cox         - 比例风险 用一个或多个解释变量预测一个事件(死亡、失败或旧病复发)发生的时间
* 时间序列     - 对误差项相关的时间序列数据建模
* 非线性       - 用一个或多个量化的解释变量预测一个量化的响应变量,不过模型是非线性的
* 非参数       - 用一个或多个量化的解释变量预测一个量化的响应变量,模型的形式源自数据形式,不事先设定
* 稳健        - 用一个或多个量化的解释变量预测一个量化的响应变量,能抵御强影响点的干扰

## symbol
| symbol        | usage   |
| --------    | :----:  |
|~  | 分隔符号,左边为响应变量,右边为解释变量。例如, 要通过x, z, w预测y,代码为 y ~ x + z + w |
|+  | 分隔预测变量 |
|:  | 表示预测变量的交互项。例如,要通过 x 、 z 及 x 与 z 的交互项预测 y ,代码为 y ~ x + z + x:z |
|*  | 表示所有可能交互项的简洁方式。代码 y~ x * z * w 可展开为 y ~ x + z + w + x:z + x:w + z:w + x:z:w |
|^  | 表示交互项达到某个次数。代码 y ~ (x + z + w)^2 可展开为 y ~ x + z + w + x:z + x:w + z:w |
|.  | 表示包含除因变量外的所有变量。例如,若一个数据框包含变量 x、y、z 和 w ,代码 y ~ . 可展开为 y ~ x + z + w |
|-  | 减号,表示从等式中移除某个变量。例如, y ~ (x + z + w)^2 – x:w 可展开为 y ~ x + z + w + x:z + z:w |
|-1  | 删除截距项。例如,表达式 y ~ x - 1 拟合 y 在 x 上的回归,并强制直线通过原点 |
|I()  | 从算术的角度来解释括号中的元素。例如, y ~ x + (z + w)^2 将展开为 y ~ x + z + w + z:w 。相反, 代码 y ~ x + I((z + w)^2 )将展开为 y ~ x + h , h 是一个由 z 和 w 的平方和创建的新变量|
|function  | 可以在表达式中用的数学函数。例如, log(y) ~ x + z + w 表示通过 x 、 z 和 w 来预测 log(y) |

## functions
* summary()      - 展示拟合模型的详细结果
* coefficients() - 列出拟合模型的模型参数(截距项和斜率)
* confint()      - 提供模型参数的置信区间(默认 95%)
* fitted()       - 列出拟合模型的预测值
* residuals()    - 列出拟合模型的残差值
* anova()        - 生成一个拟合模型的方差分析表,或者比较两个或更多拟合模型的方差分析表
* vcov()         - 列出模型参数的协方差矩阵
* AIC()          - 输出赤池信息统计量
* plot()         - 生成评价拟合模型的诊断图
* predict()      - 用拟合模型对新的数据集预测响应变量值
