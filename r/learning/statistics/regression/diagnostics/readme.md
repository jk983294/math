## 回归诊断
* 正态Q-Q图是在正态分布对应的值下,标准化残差的概率图。若满足正态假设,那么图上的点应该落在呈45度角的直线上.
* 残差图与拟合图中可以清楚地看到一个曲线关系,这暗示着你可能需要对回归模型加上一个二次项.
* 满足不变方差假设,那么在“位置尺度图”中,水平线周围的点应该随机分布.

##### “残差与杠杆图”提供了你可能关注的单个观测点的信息. 从图形可以鉴别出离群点、高杠杆值点和强影响点
* 一个观测点是离群点,表明拟合回归模型对其预测效果不佳(产生了巨大的或正或负的残差)
* 一个观测点有很高的杠杆值,表明它是一个异常的预测变量值的组合。也就是说,在预测变量空间中,它是一个离群点。因变量值不参与计算一个观测点的杠杆值。
* 一个观测点是强影响点(influential observation),表明它对模型参数的估计产生的影响过大,非常不成比例。强影响点可以通过Cook距离即Cook’s D统计量来鉴别。